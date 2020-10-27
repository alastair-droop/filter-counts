use structopt::StructOpt;
use std::path::PathBuf;
use std::process;
use signal_hook;
use std::io::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use log::*;

// Define the CLI:
#[derive(StructOpt)]
#[structopt(about="Filter HTSeq counts matrix files")]
#[structopt(after_help="Three possible filters are applied to each gene in order:
1. The gene is filtered on total read count (if -m is specified);
2. The gene is filtered on zero count (if -z is specified);
3. The gene is filtered on non-zero variance (if -i is specified)
Only genes passing all active filters are returned.")]
struct Cli {
    #[structopt(short="v", long="verbose", parse(from_occurrences), help="Provide verbose output. supply multiple times to increase verbosity")]
    verbose: usize,
    #[structopt(short="m", long="min-count", value_names=&["n"], help="Minimum total gene count")]
    min_count: Option<u32>,
    #[structopt(short="z", long="max-zerocount", value_names=&["n"], help="Maximum number of zero counts permitted in a single gene")]
    max_zero: Option<u32>,
    #[structopt(short="i", long="filter-identical", help="Filter out genes with zero variance (i.e. with all values identical)")]
    filter_identical: bool,
    #[structopt(short="o", long="output-metacounts", help="Extract metacounts (starting with double underscores) to file")]
    metacount_path: Option<PathBuf>,
    #[structopt(parse(from_os_str), help="Input counts file")]
    path: PathBuf,
}

fn main() {
    // Capture the commandline arguments:
    let args = Cli::from_args();
    
    // Register the SIGPIPE actions:
    let _signal = unsafe { signal_hook::register(signal_hook::SIGPIPE, || process::exit(128 + signal_hook::SIGPIPE)) };
        
    // Build the log:
    stderrlog::new().module(module_path!()).verbosity(args.verbose).init().expect("Failed to initialise logging");

    // Attempt to open the input file:
    info!("{}", format!("reading counts from {}", args.path.to_str().expect("Failed to extract input path")));
    let input_file = File::open(&args.path).expect("Failed to open input file");
    let input_buffer = BufReader::new(input_file);
    
    // If necessary, attempt to open the metacount output file:
    let mut metacount_buffer: Option<BufWriter<std::fs::File>> = None;
    match args.metacount_path {
        Some(p) => {
            info!("{}", format!("saving metacounts to {}", p.to_str().expect("Failed to extract metacount path")));
            let metacount_file = File::create(&p).expect("Failed to open metacount file");
            metacount_buffer = Some(BufWriter::new(metacount_file));
        },
        _ => (),
    }
    
    // Get the line iterator:
    let mut line_iter = input_buffer.lines().map(|l| l.expect("Failed to read line from input file"));

    // Grab the header:
    let header = line_iter.next();
    match header {
        Some(h) => {
            println!("{}", h);
            match metacount_buffer {
                Some(mut b) => {
                    writeln!(b, "{}", h);
                },
                _ => (),
            }
        },
        None => {
            eprintln!("No header line detected");
            process::exit(1);
        },
    };
    
    // Record the total & filtered genes:
    let mut total_genes: u32 = 0;
    let mut passed_genes: u32 = 0;
    let mut metagenes: u32 = 0;
    
    // Iterate over the remaining lines:
    for line in line_iter {
        // Pull out the line, and split it into parts:
        let line_trimmed = line.trim();
        let line_data: Vec<&str> = line_trimmed.split('\t').collect();

        // Check if this is a metacount gene; if so write it to the output file and move on:
        if line_data[0].starts_with("__") {
            metagenes += 1;
            debug!("{}", format!("metacount {} detected", line_data[0]));
            match metacount_buffer {
                Some(mut b) => {
                    writeln!(b, "{}", line_trimmed.trim_start_matches("_"));
                    continue;
                },
                _ => (),
            }
        } else {
            // Record the gene:
            total_genes += 1;                
        }

        let mut total: u32 = 0;
        let mut all_equal: bool = true;
        let mut n_zero = 0;
        let mut last_value: u32 = line_data[1].parse::<u32>().expect("Failed to parse counts line");
        for i in line_data[1..].iter() {
            let x = i.parse::<u32>().expect("Failed to parse counts line");
            total += x;
            if x == 0 {
                n_zero += 1;
            }
            if x != last_value {
                all_equal = false;
            }
            last_value = x;
        }
        // Filter line by minimum count:
        match args.min_count {
            Some(i) if total < i => {
                debug!("{}", format!("gene {} filtered (total count {} < {})", line_data[0], total, i));
                continue
            },
            _ => (),
        }
        // Filter on maximum zero count:
        match args.max_zero {
            Some(i) if n_zero > i => {
                debug!("{}", format!("gene {} filtered (zero count {} > {})", line_data[0], n_zero, i));
                continue
            },
            _ => (),
        }
        // Filter on non-zero variance:
        if all_equal && args.filter_identical {
            debug!("{}", format!("gene {} filtered (zero variance)", line_data[0]));
            continue
        };
        // Record the gene passed:
        passed_genes += 1;
        // Output the (non-filtered) line:
        println!("{}", line_trimmed);
    }
    
    // Record the results:
    info!("{}", format!("{} / {} genes passed filter", passed_genes, total_genes));
    info!("{}", format!("{} metagenes detected", metagenes));
    
}

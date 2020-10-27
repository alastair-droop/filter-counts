use structopt::StructOpt;
use std::path::PathBuf;
use std::process;
use signal_hook;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;
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

    // Get the line iterator:
    let mut line_iter = input_buffer.lines().map(|l| l.expect("Failed to read line from input file"));

    // Grab the header:
    let header = line_iter.next();
    match header {
        Some(h) => {
            println!("{}", h);
        },
        None => {
            eprintln!("No header line detected");
            process::exit(1);
        },
    };
    
    // Record the total & filtered genes:
    let mut total_genes: u32 = 0;
    let mut passed_genes: u32 = 0;
    
    // Iterate over the remaining lines:
    for line in line_iter {
        // Record the gene:
        total_genes += 1;
        // Pull out the line, and split it into parts:
        let line_trimmed = line.trim();
        let line_data: Vec<&str> = line_trimmed.split('\t').collect();
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
}

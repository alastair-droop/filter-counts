use structopt::StructOpt;
use std::path::PathBuf;
use std::process;
use signal_hook;
use std::io::prelude::*;
use std::fs::File;
use std::io::{Error, ErrorKind};
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
    min_count: Option<usize>,
    #[structopt(short="z", long="max-zerocount", value_names=&["n"], help="Maximum number of zero counts permitted in a single gene")]
    max_zero: Option<usize>,
    #[structopt(short="i", long="filter-identical", help="Filter out genes with zero variance (i.e. with all values identical)")]
    filter_identical: bool,
    #[structopt(short="o", long="output-metacounts", help="Extract metacounts (starting with double underscores) to file")]
    metacount_path: Option<PathBuf>,
    #[structopt(parse(from_os_str), help="Input counts file")]
    path: PathBuf,
}

fn main() -> Result<(), Error> {
    // Capture the commandline arguments:
    let args = Cli::from_args();

    // Register the SIGPIPE actions:
    let _signal = unsafe { signal_hook::register(signal_hook::SIGPIPE, || process::exit(128 + signal_hook::SIGPIPE)) };

    // Build the log:
    if stderrlog::new().module(module_path!()).verbosity(args.verbose).init().is_err() {
        return Err(Error::new(ErrorKind::Other, "failed to initialise logger"));
    }

    // Attempt to open the input file:
    let input_file = File::open(&args.path)?;
    if let Some(file_name) = args.path.to_str() {
        info!("{}", format!("reading counts from {}", file_name));
    }
    let input_buffer = BufReader::new(input_file);
    
    // If necessary, attempt to open the metacount output file:
    let mut metacount_buffer = match args.metacount_path {
        Some(p) => {
            let f = File::create(&p)?;
            if let Some(f_name) = p.to_str() {
                info!("{}", format!("saving metacounts to {}", f_name));
            }
            Some(BufWriter::new(f))
        },
        None => None,
    };
    
    // Get the line iterator:
    let mut line_iter = input_buffer.lines();
    
    // Read the file header:
    let file_header_result = match line_iter.next() {
        Some(h) => h,
        None => Err(Error::new(ErrorKind::UnexpectedEof, "failed to read input file header")),
    };
    let file_header = match file_header_result {
        Ok(h) => h,
        Err(_) => return Err(Error::new(ErrorKind::InvalidData, "failed to parse file header")),
    };
    
    // Write out the file header:
    println!("{}", file_header);
    if let Some(b) = &mut metacount_buffer {
        writeln!(b, "{}", file_header)?;
    }
    
    // Record the total & filtered genes:
    let mut total_genes: usize = 0;
    let mut passed_genes: usize = 0;
    let mut metagenes: usize = 0;

    // Iterate over the remaining lines:
    for line_res in line_iter {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => return Err(Error::new(ErrorKind::InvalidData, "failed to parse input file")),
        };
        let line_trimmed = line.trim();
        let line_data: Vec<_> = line_trimmed.split('\t').collect();
        let gene = &line_data[0];
        // Check if this is a metagene:
        if gene.starts_with("__") {
            metagenes += 1;
            if let Some(b) = &mut metacount_buffer {
                writeln!(b, "{}", line_trimmed.trim_start_matches('_'))?;
                debug!("{}", format!("metacount \"{}\" written to separate output", gene));
            } else {
                println!("{}", line_trimmed);
                debug!("{}", format!("metacount \"{}\" included in main output", gene));
            }
            continue;
        }
        let counts: Result<Vec<_>, _> = line_data.iter().skip(1).map(|s| s.parse::<usize>()).collect();
        match counts {
            Ok(c) => {
                total_genes += 1;
                // Filter on sum:
                if let Some(min_count) = args.min_count {
                    let line_sum: usize = c.iter().sum();
                    if line_sum < min_count{
                        debug!("{}", format!("gene {} filtered (total count {} < {})", gene, line_sum, min_count));
                        continue;
                    }
                }
                // Filter on zero count:
                if let Some(max_zcount) = args.max_zero {
                    let line_zcount: usize = c.iter().filter(|i| **i == 0).count();
                    if line_zcount > max_zcount {
                        debug!("{}", format!("gene {} filtered (zero count {} > {})", gene, line_zcount, max_zcount));
                        continue;
                    }
                }
                // Filter on zero variance:
                if args.filter_identical && c.iter().all(|i| *i == c[0]) {
                    debug!("{}", format!("gene {} filtered (zero variance)", gene));
                    continue;
                }
                // Gene passed filtering:
                passed_genes += 1;
                debug!("{}", format!("gene {} passed filtering", gene));
                println!("{}", line_trimmed);
            },
            Err(_) => {
                warn!("{}", format!("failed to convert counts from line {}", line.trim()));
                continue;
            },
        }
    }
    // Record the results:
    info!("{}", format!("{} / {} genes passed filter", passed_genes, total_genes));
    info!("{}", format!("{} metagenes detected", metagenes));
    Ok(())
}

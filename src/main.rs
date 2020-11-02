use shellexpand;
use structopt::StructOpt;
use std::path::PathBuf;
use std::process;
use signal_hook;
use std::io::prelude::*;
use std::fs::File;
use std::io::{stdout, Error, ErrorKind};
use std::io::{BufReader, BufWriter};
use log::*;

// Define a struct to hold sample metadata:
#[derive(Debug)]
struct Sample {
    name: String,
    metacounts: Vec<u64>,
    total_count: u64,
    passed_count: u64,
    total_expressed: u64,
    passed_expressed: u64,
}

// Define a struct to record how we're outputting metacounts:
struct MetacountDestination {
    handle: Box<dyn Write>,
    is_stdout: bool,
    prefix: String,
}

// A function to expand a PathBuf to a full string:
fn expand_path(path: &PathBuf) -> Option<String> {
    let path_str = match path.to_str() {
        Some(p) => p,
        None => return None,
    };
    let path_exp = match shellexpand::full(path_str) {
        Ok(p) => p,
        Err(_) => return None,
    };
    Some(String::from(path_exp))
}

// Define the CLI:
#[derive(StructOpt)]
#[structopt(about="Filter HTSeq counts matrix files")]
#[structopt(after_help="The following filters are applied to each gene:
* The gene is filtered on total read count (if -m is specified);
* The gene is filtered on zero count (if -z is specified);
* The gene is filtered on non-zero variance (if -i is specified)")]
struct Cli {
    #[structopt(short="v", long="verbose", parse(from_occurrences), help="Provide verbose output. supply multiple times to increase verbosity")]
    verbose: usize,
    #[structopt(short="m", long="min-count", value_names=&["n"], help="Minimum total gene count")]
    min_count: Option<u64>,
    #[structopt(short="e", long="min-expressed", value_names=&["n"], help="Minimum number of expressed samples")]
    min_expressed: Option<u64>,
    #[structopt(short="i", long="filter-identical", help="Filter out genes with zero variance (i.e. with all values identical)")]
    filter_identical: bool,
    #[structopt(short="x", long="expression", value_names=&["e"], default_value="1", help="Minimum expression count")]
    expression_threshold: u64,
    #[structopt(parse(from_os_str), short="o", long="metacount-file", value_names=&["path"], help="Extract metacounts (starting with double underscores) to file")]
    metacount_path: Option<PathBuf>,
    #[structopt(short="s", long="summary", help="Include sample summary metacounts")]
    summary_metacounts: bool,
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
    let input_filename = match expand_path(&args.path) {
        Some(f) => f,
        None => return Err(Error::new(ErrorKind::NotFound, "input file not found")),
    };
    let input_file = File::open(&input_filename)?;
    let input_buffer = BufReader::new(input_file);
    info!("{}", format!("reading counts from {}", input_filename));    
    
    // Sort out the metacount destination:
    let mut metacount_dest = match args.metacount_path {
        Some(p) => {
            let metacount_filename = match expand_path(&p) {
                Some(f) => f,
                None => return Err(Error::new(ErrorKind::NotFound, "metacount output file not found")),
            };
            let f = File::create(&metacount_filename)?;
            info!("{}", format!("writing metacounts to {}", metacount_filename));
            MetacountDestination{
                handle: Box::new(BufWriter::new(f)) as Box<dyn Write>,
                is_stdout: false,
                prefix: String::from(""),
            }
        },
        None => {
            info!("writing metacounts to stdout");
            MetacountDestination{
                handle: Box::new(stdout()) as Box<dyn Write>,
                is_stdout: true,
                prefix: String::from("__"),
            }
        }
    };
    
    // Assign a Vec to capture the metacount names:
    let mut metacount_names: Vec<String> = Vec::with_capacity(5);

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

    // Initialise the sample metadata structs:
    let mut samples: Vec<Sample> = Vec::new();
    for sample in file_header.trim().split('\t').skip(1) {
        samples.push(Sample{
            name: String::from(sample),
            metacounts: Vec::with_capacity(5),
            total_count: 0,
            passed_count: 0,
            total_expressed: 0,
            passed_expressed: 0,
        })
    }

    // Record the total & filtered genes:
    let mut total_genes: u64 = 0;
    let mut passed_genes: u64 = 0;

    // Iterate over the remaining lines:
    for line_res in line_iter {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => return Err(Error::new(ErrorKind::InvalidData, "failed to parse input file")),
        };
        let line_trimmed = line.trim();
        let line_data: Vec<_> = line_trimmed.split('\t').collect();
        let gene = &line_data[0];

        // Extract the counts:
        let counts: Vec<_> = match line_data.iter().skip(1).map(|s| s.parse::<u64>()).collect() {
            Ok(c) => c,
            Err(_) => {
                warn!("{}", format!("failed to convert counts from line {}", line.trim()));
                continue;
            }
        };

        // Check if this is a metagene:
        if gene.starts_with("__") {
            metacount_names.push(String::from(*gene));
            for m in counts.iter().enumerate(){
                samples[m.0].metacounts.push(*m.1);
            }
            continue;
        }

        total_genes += 1;

        // Calculate the gene stats:
        let mut gene_total: u64 = 0;
        let mut gene_nexpressed: u64 = 0;
        let mut gene_filtered = false;
        for i in counts.iter().enumerate() {
            gene_total += i.1;
            samples[i.0].total_count += i.1;
            if *i.1 >= args.expression_threshold {
                gene_nexpressed += 1;
                samples[i.0].total_expressed += 1;
            }
        }

        // Filter by minimum count:
        match args.min_count {
            Some(min_count) if gene_total < min_count => {
                debug!("{}", format!("gene {} failed filtering (total count {} < {})", gene, gene_total, min_count));
                gene_filtered = true;
            },
            _ => (),
        };

        // Filter on zero count:
        match args.min_expressed {
            Some(min_expressed) if gene_nexpressed < min_expressed => {
                debug!("{}", format!("gene {} failed filtering (expressed count {} < {})", gene, gene_nexpressed, min_expressed));
                gene_filtered = true;
            },
            _ => (),
        }

        // Filter on zero variance:
        if args.filter_identical && counts.iter().all(|i| *i == counts[0]) {
            debug!("{}", format!("gene {} failed filtering (zero variance)", gene));
            gene_filtered = true;
        }

        if !gene_filtered {
            // Gene passed filtering:
            passed_genes += 1;
            trace!("{}", format!("gene {} passed filtering", gene));
            println!("{}", line_trimmed);
            for i in counts.iter().enumerate() {
                samples[i.0].passed_count += i.1;
                if *i.1 >= args.expression_threshold {
                    samples[i.0].passed_expressed += 1;
                }
            }
        }
    }

    // Process and write the metacount data:
    if !metacount_dest.is_stdout {
        let mut header: Vec<String> = vec!("feature".to_string());
        header.extend(samples.iter().map(|s|s.name.to_string()));
        writeln!(metacount_dest.handle, "{}", header.join("\t"))?;        
    }
    if args.summary_metacounts {
        writeln!(metacount_dest.handle, "{}total_count\t{}", metacount_dest.prefix, samples.iter().map(|s|s.total_count.to_string()).collect::<Vec<String>>().join("\t"))?;
        writeln!(metacount_dest.handle, "{}passed_count\t{}", metacount_dest.prefix, samples.iter().map(|s|s.passed_count.to_string()).collect::<Vec<String>>().join("\t"))?;
        writeln!(metacount_dest.handle, "{}total_expressed\t{}", metacount_dest.prefix, samples.iter().map(|s|s.total_expressed.to_string()).collect::<Vec<String>>().join("\t"))?;
        writeln!(metacount_dest.handle, "{}passed_expressed\t{}", metacount_dest.prefix, samples.iter().map(|s|s.passed_expressed.to_string()).collect::<Vec<String>>().join("\t"))?;
    }
    for m in metacount_names.iter().enumerate() {
        let counts = samples.iter().map(|s|s.metacounts[m.0].to_string()).collect::<Vec<String>>().join("\t");
        if metacount_dest.is_stdout {
            writeln!(metacount_dest.handle, "{}\t{}", m.1, counts)?;
        } else {
            writeln!(metacount_dest.handle, "{}\t{}", m.1.trim_start_matches('_'), counts)?;
        }
    }

    // Record the results:
    info!("{}", format!("{} / {} genes passed filter", passed_genes, total_genes));
    info!("{}", format!("{} metagenes detected", metacount_names.len()));
    Ok(())
}

mod constants;
mod dna_helpers;
mod helpers;
mod output;
mod train;
mod viterbi;

use bio::io::fasta;
use clap::{App, Arg};
use dna_helpers::get_prob_from_cg;
use helpers::create_file_if_not_exists;
use output::print_prediction;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self};
use std::path::Path;
use train::{Train, HMM};
use viterbi::{viterbi};
use std::sync::{Arc, Mutex};


/**
 * main.rs
 * =======
 * This file contains the main method used for calling the viterbi function from the command line
 * and processing the output.
 */

fn main() {
    // We use clap to process command line arguments
    let matches = App::new("fgsrs")
        .version("1.0")
        .author("Laurens Debackere <Laurens.Debackere@UGent.be>")
        .about("A reimplementation of the original FragGeneScan project (see https://omics.informatics.indiana.edu/FragGeneScan/) in Rust with an improved command-line interface, better performance and code quality.")
        .arg(
            Arg::with_name("whole-genome")
                .short("w")
                .long("whole-genome")
                .value_name("WHOLE_GENOME")
                .help("(REQUIRED) Does the input contain sequence reads (WHOLE_GENOME = 0) or full genome sequences (WHOLE_GENOME = 1).")
                .takes_value(true)
                .required(true)
        )
        .arg(
            Arg::with_name("threads")
                .short("p")
                .long("threads")
                .value_name("NUM_THREADS")
                .help("How many threads should the program use.")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("train")
                .short("t")
                .long("train")
                .value_name("TRAIN_PATH")
                .help("(REQUIRED) Path to a directory containing the train model files to be used.")
                .takes_value(true)
                .required(true)
        )
        .arg(
            Arg::with_name("output")
                .short("d")
                .long("output")
                .value_name("DNA_OUTPUT_FILE")
                .help("(OPTIONAL) Specifies a file path where the DNA-FASTA-file is written to.")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("metadata")
                .short("e")
                .long("metadata")
                .value_name("_OUTPUT_FILE")
                .help("(OPTIONAL) Specifies a file path where the Metadata-file is written to.")
                .takes_value(true)
        )
        .get_matches();

    /*
    Convert the specified training path into a Train struct
    */
    let train_path = matches.value_of("train").unwrap();
    let train;
    let hmm;
    if Path::new(train_path).exists() {
        train = Train::from_file();
        hmm = HMM::from_file(train_path);
    } else {
        println!("ERROR: Something went wrong while accessing the specified training directory.");
        return;
    }

    /*
    Process the -w parameter
    */
    let wholegenome = matches.value_of("whole-genome").unwrap();
    let wholegenome = match wholegenome {
        "1" => true,
        "0" => false,
        _ => {
            println!("ERROR: Invalid value specified for -w.");
            return;
        }
    };

    /*
    Process the -p parameter
    */
    if matches.is_present("threads") {
        let threadnum: usize = matches
            .value_of("threads")
            .unwrap()
            .parse()
            .expect("ERROR: The parameter -p should have a numeric value.");
        if threadnum < 1 {
            println!("ERROR: Invalid number of threads");
            return;
        }
        // Configure threadpool based on specifed number of threads
        rayon::ThreadPoolBuilder::new()
            .num_threads(threadnum)
            .build_global()
            .unwrap();
    }

    /*
     * Next we should read the fasta sequences from STDIN, we will only process those sequences longer than 70 bp's
     */
    let records: Vec<(String, String)> = fasta::Reader::new(io::stdin())
        .records()
        .map(|result| {
            let record = result.unwrap();
            // obtain sequence
            let seq = std::str::from_utf8(record.seq()).unwrap().to_string();
            let id = String::from(record.id());
            (id, seq)
        })
        .filter(|(_, seq)| seq.len() > 70)
        .collect();

    
    /*
     * Process -e parameter to get output metadata file
     */
    let mut metadata_output: Arc<Mutex<Option<File>>> = Arc::new(Mutex::new(None));
    if matches.is_present("metadata") {
        metadata_output = Arc::new(Mutex::new(Some(create_file_if_not_exists(
            matches.value_of("metadata").unwrap(),
        ))));
    }

    /*
     * Process -d parameter to get output dna file
     */
    let mut dna_output: Arc<Mutex<Option<File>>> = Arc::new(Mutex::new(None));
    if matches.is_present("output") {
        dna_output = Arc::new(Mutex::new(Some(create_file_if_not_exists(
            matches.value_of("output").unwrap(),
        ))));
    }

    /*
     * Now we use the Rayon Parallel Iterator to perform the viterbi algorithm on each of the sequences,
     * this will return a Vec of predictions
     */
    records
        .into_par_iter()
        .for_each(|(header, sequence)| {
            let mut local_hmm = hmm.clone();
            let cg = get_prob_from_cg(&mut local_hmm, &train, &sequence);
            let pred = viterbi(&local_hmm, &train, &sequence, wholegenome, cg, &header);
            print_prediction(pred, &metadata_output, &dna_output);
        });
}

mod constants;
mod train;
mod viterbi;
mod helpers;
mod output;

use bio::io::fasta;
use clap::{App, Arg};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self};
use std::path::Path;
use train::{get_prob_from_cg, Train, HMM};
use viterbi::{viterbi, Prediction};
use output::{print_prediction};
use helpers::{create_file_if_not_exists};

fn main() {
    let matches = App::new("FragGeneScan.Rs")
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
        rayon::ThreadPoolBuilder::new()
            .num_threads(threadnum)
            .build_global()
            .unwrap();
    }

    /*
     * Next we should read the fasta sequences from STDIN, start our threads and do some work
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

    let predictions: Vec<Prediction> = records
        .into_par_iter()
        .map(|(header, sequence)| {
            let mut local_hmm = hmm.clone();
            let cg = get_prob_from_cg(&mut local_hmm, &train, &sequence);
            viterbi(&local_hmm, &train, &sequence, wholegenome, cg, &header)
        })
        .collect();

    /*
     * Join all of the work of these threads
     */

    let mut metadata_output: Option<File> = None;
    if matches.is_present("metadata") {
        metadata_output = Some(create_file_if_not_exists(
            matches.value_of("metadata").unwrap(),
        ));
    }

    let mut dna_output: Option<File> = None;
    if matches.is_present("output") {
        dna_output = Some(create_file_if_not_exists(
            matches.value_of("output").unwrap(),
        ));
    }

    predictions.into_iter().for_each(|prediction| {
        print_prediction(prediction, &mut metadata_output, &mut dna_output);
        }
    );
}




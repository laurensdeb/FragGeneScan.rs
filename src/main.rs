use bio::io::fasta;
use std::io::{self, Write};
use std::path::Path;
mod types;
use self::types::train::{get_prob_from_cg, Train, HMM};
use self::types::viterbi::{viterbi, Output};
use rayon::prelude::*;
use clap::{App, Arg};
use std::fs::File;



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
                .default_value("1")
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
    let wholegenome = match wholegenome { "1" => true, "0" => false, _ => {
        println!("ERROR: Invalid value specified for -w.");
        return;
    }};

    /*
    Process the -p parameter
    */
    let threadnum: usize = matches.value_of("threads").unwrap().parse().expect("ERROR: The parameter -t should have a numeric value.");
    if threadnum < 1 {
        println!("ERROR: Invalid number of threads");
        return;
    }

    /*
     * Next we should read the fasta sequences from STDIN, start our threads and do some work
     */
    let format = false; // TODO: this parameter was undocumented up until now...
    let reader = fasta::Reader::new(io::stdin());
    let mut records: Vec<(String, String)> = Vec::new();


    for result in reader.records(){
        // obtain record or fail with error
        let record = result.unwrap();
        // obtain sequence
        let seq = record.seq();
        let desc = String::from(record.id());
        let sequence = std::str::from_utf8(seq).unwrap().to_string();
        records.push((desc, sequence))
    }

    let outputs: Vec<Output> = records.into_par_iter()
    .map(|(header, sequence)| {
        let mut local_hmm = hmm.clone();
        let cg = get_prob_from_cg(&mut local_hmm, &train, &sequence);
        viterbi(
            &local_hmm,
            &train,
            &sequence,
            wholegenome,
            cg,
            format,
            &header,
        )
    }).collect();

     /*
     * Join all of the work of these threads
     */ 
    let out = io::stdout();
    let mut handle = out.lock();
    for output in outputs{
        write!(handle, "{}", output.aa);
    }
    drop(handle);

     /*
     * Write output to optional files and stdout
    */

    /*
    // TODO: Process command line arguments
    // TODO: fill train and hmm structs using cmd arguments
    let reader = fasta::Reader::new(io::stdin());
    let _train = Train::from_file();
    let _hmm = HMM::from_file(&String::from("train/illumina_1"));
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        // obtain sequence
        let seq = record.seq();
        let string = std::str::from_utf8(seq).unwrap().to_string();
        println!("{}", string);
    }
    */
}
/*
pub struct ThreadData {
    out_file: File,
    aa_file: File,
    dna_file: File,
    head: String,
    sequence: String,
    wholegenome: bool,
    cg: usize,
    format: bool,
    hmm: HMM,
    train: Train,
}

fn thread_func(data: &mut ThreadData) {
    data.cg = get_prob_from_cg(&mut data.hmm, &data.train, &data.sequence); //cg - 26 Ye April 16, 2016
    if data.sequence.chars().count() > 70 {
        viterbi(
            &data.hmm,
            &data.train,
            &data.sequence,
            data.wholegenome,
            data.cg,
            data.format,
            &data.head,
            &mut data.out_file,
            &mut data.aa_file,
            &mut data.dna_file,
        )
    }
}
*/
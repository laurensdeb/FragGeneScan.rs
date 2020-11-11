use super::helpers::write_data;
use super::viterbi::Prediction;
use std::fs::File;
use std::sync::{Arc, Mutex};


/**
 * output.rs
 * =========
 * This file contains methods responsible for processing the Prediction structs resulting from a run of viterbi::viterbi,
 * it will take care of writing the output to the proper channel (stdout / (optional) metadata files).
 */

/**
 * This method will print a single Prediction to stdout as AA sequence and optionally output metadata and dna metadata.
 */
pub fn print_prediction(
    prediction: Prediction,
    metadata_output: &Arc<Mutex<Option<File>>>,
    dna_output: &Arc<Mutex<Option<File>>>,
) {
    let metadata_option = &mut *(metadata_output.lock().unwrap());
    let dna_option = &mut *(dna_output.lock().unwrap());

    
    // Should we output to the metadata file
    if metadata_option.is_some() {
        write_data(
            metadata_option.as_mut().unwrap(),
            format!(">{}\n", prediction.head),
        );
    }
    for out in prediction.outs {
        print_aa(
            &prediction.head,
            out.dna_start_t,
            out.dna_end_t,
            out.forward,
            &out.protein,
        );
        // Should we output to the metadata file
        if metadata_option.is_some() {
            print_metadata(
                metadata_option.as_mut().unwrap(),
                out.dna_start_t,
                out.dna_end_t,
                out.frame,
                out.forward,
                out.final_score,
                out.insert,
                out.delete,
            );
        }
        // Should we output to the dna metadata file
        if dna_option.is_some() {
            print_dna_metadata(
                dna_option.as_mut().unwrap(),
                &prediction.head,
                out.dna_start_t,
                out.dna_end_t,
                out.forward,
                &out.dna,
            );
        }
    }
}

/**
 * Helper method to print amino acids to stdout in correct format
 */
fn print_aa(head: &String, start_t: usize, end_t: usize, forward: bool, protein: &String) {
    println!(
        ">{}_{}_{}_{}",
        head,
        start_t,
        end_t,
        forward_to_chr(forward)
    );
    println!("{}", protein);
}

/**
 * Helper method to write dna output metadata to the specified output file
 */
fn print_dna_metadata(
    dna_output: &mut File,
    head: &String,
    start_t: usize,
    end_t: usize,
    forward: bool,
    dna: &String,
) {
    write_data(
        dna_output,
        format!(
            ">{}_{}_{}_{}\n",
            head,
            start_t,
            end_t,
            forward_to_chr(forward)
        ),
    );
    write_data(dna_output, format!("{}\n", dna));
}

/**
 * Helper method to write metadata to the specified output file
 */
fn print_metadata(
    out: &mut File,
    start_t: usize,
    end_t: usize,
    frame: usize,
    forward: bool,
    final_score: f64,
    insert: Vec<usize>,
    delete: Vec<usize>,
) {
    write_data(
        out,
        format!("{}\t{}\t{}\t{}\t{}\t", start_t, end_t, forward_to_chr(forward), frame, final_score),
    );
    write_data(out, String::from("I:"));
    for i in insert {
        write_data(out, format!("{},", i));
    }
    write_data(out, String::from("\tD:"));
    for d in delete {
        write_data(out, format!("{},", d));
    }
    write_data(out, String::from("\n"));
}

/**
 * Helper method to convert the strand (forward/reverse) to a +/- character respectively.
 */
fn forward_to_chr(forward: bool) -> char {
    match forward {
        true => '+',
        false => '-',
    }
}

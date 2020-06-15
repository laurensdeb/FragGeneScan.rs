use super::viterbi::{Prediction};
use super::helpers::{write_data};
use std::fs::File;


pub fn print_prediction(prediction: Prediction, metadata_output: &mut Option<File>, dna_output: &mut Option<File>) {
    if metadata_output.is_some() {
        write_data(
            metadata_output.as_mut().unwrap(),
            format!("{}\n", prediction.head),
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
        if metadata_output.is_some() {
            print_metadata(
                metadata_output.as_mut().unwrap(),
                out.dna_start_t,
                out.dna_end_t,
                out.frame,
                out.final_score,
                out.insert,
                out.delete,
            );
        }
        if dna_output.is_some() {
            print_dna_metadata(
                dna_output.as_mut().unwrap(),
                &prediction.head,
                out.dna_start_t,
                out.dna_end_t,
                out.forward,
                &out.dna,
            );
        }
    }
}


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

fn print_metadata(
    out: &mut File,
    start_t: usize,
    end_t: usize,
    frame: usize,
    final_score: f64,
    insert: Vec<usize>,
    delete: Vec<usize>,
) {
    write_data(
        out,
        format!("{}\t{}\t+\t{}\t{}\t", start_t, end_t, frame, final_score),
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

fn forward_to_chr(forward: bool) -> char {
    match forward {
        true => '+',
        false => '-',
    }
}
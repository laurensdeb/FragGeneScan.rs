use super::train::{Train, HMM};
use rayon::prelude::*;

/**
 * dna_helpers.rs
 * ==============
 * This file contains a number of helper methods related to processing DNA/AA values
 * from string to integer and vice versa as well as some methods related to choosing the
 * correct training values based on CG-score.
 */

const CODON: &'static [char] = &['A', 'C', 'G', 'T', 'N'];

const CODON_CODE: &'static [char] = &[
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q',
    'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A',
    'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*',
    'C', 'W', 'C', 'L', 'F', 'L', 'F', 'X',
];

const ANTI_CODON_CODE: &'static [char] = &[
    'F', 'V', 'L', 'I', 'C', 'G', 'R', 'S', 'S', 'A', 'P', 'T', 'Y', 'D', 'H', 'N', 'L', 'V', 'L',
    'M', 'W', 'G', 'R', 'R', 'S', 'A', 'P', 'T', '*', 'E', 'Q', 'K', 'F', 'V', 'L', 'I', 'C', 'G',
    'R', 'S', 'S', 'A', 'P', 'T', 'Y', 'D', 'H', 'N', 'L', 'V', 'L', 'I', '*', 'G', 'R', 'R', 'S',
    'A', 'P', 'T', '*', 'E', 'Q', 'K', 'X',
];

/*
Code for loading the relevant information into the HMM struct for this CG value of the sequence
*/
pub fn get_prob_from_cg(hmm: &mut HMM, train: &Train, seq: &String) -> usize {
    //change from void to int, Ye, April 18, 2016
    let mut cg_count: usize = seq.chars().collect::<Vec<char>>().into_par_iter().map(|c|  match c {
            'c' | 'C' | 'g' | 'G' => 1,
            _ => 0,
        }).sum();
    
    let len_seq = seq.len();

    let cg_count_i = ((((cg_count as f64 * 1.0) / len_seq as f64) * 100.0).floor() as i32) - 26;
    if cg_count_i < 0 {
        cg_count = 0;
    } else if cg_count_i > 43 {
        cg_count = 43;
    } else {
        cg_count = cg_count_i as usize;
    }
    hmm.e_m = train.trans[cg_count].clone();
    hmm.e_m_1 = train.rtrans[cg_count].clone();
    hmm.tr_r_r = train.noncoding[cg_count].clone();
    hmm.tr_s = train.start[cg_count].clone();
    hmm.tr_e = train.stop[cg_count].clone();
    hmm.tr_s_1 = train.start1[cg_count].clone();
    hmm.tr_e_1 = train.stop1[cg_count].clone();
    hmm.s_dist = train.s_dist[cg_count].clone();
    hmm.e_dist = train.e_dist[cg_count].clone();
    hmm.s1_dist = train.s1_dist[cg_count].clone();
    hmm.e1_dist = train.e1_dist[cg_count].clone();
    cg_count
}

/**
 * Converts a nucleotide character to an integer
 */
pub fn nt2int(nt: char) -> usize {
    match nt {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4,
    }
}

fn nt2int_rc(nt: &char) -> usize {
    match nt {
        'A' | 'a' => 3,
        'C' | 'c' => 2,
        'G' | 'g' => 1,
        'T' | 't' => 0,
        _ => 4,
    }
}

/**
 * Converts a transition specifier string to an integer
 */
pub fn tr2int(tr: &String) -> usize {
    match tr.as_str() {
        "MM" => 0,
        "MI" => 1,
        "MD" => 2,
        "II" => 3,
        "IM" => 4,
        "DD" => 5,
        "DM" => 6,
        "GE" => 7,
        "GG" => 8,
        "ER" => 9,
        "RS" => 10,
        "RR" => 11,
        "ES" => 12,
        "ES1" => 13,
        _ => panic!(format!("Unknown transition: {}", tr)),
    }
}

/**
 * Converts a trinucleotide to an integer
 */
pub fn trinucleotide(a: &char, b: &char, c: &char) -> usize {
    let mut result = match a {
        'A' | 'a' => 0,
        'C' | 'c' => 16,
        'G' | 'g' => 32,
        'T' | 't' => 48,
        _ => 0,
    };
    result += match b {
        'A' | 'a' => 0,
        'C' | 'c' => 4,
        'G' | 'g' => 8,
        'T' | 't' => 12,
        _ => 0,
    };
    result += match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 0,
    };
    result
}

fn trinucleotide_pep(a: &char, b: &char, c: &char) -> usize {
    let mut result = match a {
        'A' | 'a' => 0,
        'C' | 'c' => 16,
        'G' | 'g' => 32,
        'T' | 't' => 48,
        _ => 64,
    };
    if result < 64 {
        result += match b {
            'A' | 'a' => 0,
            'C' | 'c' => 4,
            'G' | 'g' => 8,
            'T' | 't' => 12,
            _ => {
                result = 64;
                0
            }
        };
    }
    if result < 64 {
        result += match c {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => {
                result = 64;
                0
            }
        };
    }
    result
}

/**
 * Get reverse coding DNA
 */
pub fn get_rc_dna(dna: &Vec<char>) -> Vec<char> {
    dna.par_iter().rev().map(|c| CODON[nt2int_rc(&*c)]).collect()
}

/**
 * Get protein sequence from dna sequence
 */
pub fn get_protein(dna: &Vec<char>, strand: bool, wholegenome: bool) -> Vec<char> {
    let len = dna.len();
    let mut protein = vec!['\0'; len / 3];
    if strand {
        protein.par_iter_mut().enumerate().for_each(|(i, p)| {
            let i = i * 3;
            *p = CODON_CODE[trinucleotide_pep(&dna[i], &dna[i + 1], &dna[i + 2])]
        })
    } else {
        for i in (0..len).step_by(3) {
            if (len - i) / 3 > 0 {
                protein[(len - i) / 3 - 1] =
                    ANTI_CODON_CODE[trinucleotide_pep(&dna[i], &dna[i + 1], &dna[i + 2])];
            }
        }
    }
    if protein[len / 3 - 1] == '*' {
        protein.pop();
    }
    if wholegenome {
        return protein; //short reads, skip
    }
    if strand {
        let s = trinucleotide_pep(&dna[0], &dna[1], &dna[2]);
        // trinucleotide_pep(&'G', &'T', &'G') == 46
        // trinucleotide_pep(&'T', &'T', &'G') == 62
        if s == 46 || s == 62 {
            protein[0] = 'M';
        }
    } else {
        let s = trinucleotide_pep(&dna[len - 3], &dna[len - 2], &dna[len - 1]);
        // trinucleotide_pep(&'C', &'A', &'C') == 17
        // trinucleotide_pep(&'C', &'A', &'A') == 16
        if s == 17 || s == 16 {
            protein[0] = 'M';
        }
    }
    protein
}

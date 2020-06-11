use std::fs::File;
use std::io::{BufRead, BufReader};
use whiteread::parse_string;

const MFILENAME: &'static str = "train/gene";
const M1FILENAME: &'static str = "train/rgene";
const NFILENAME: &'static str = "train/noncoding";
const SFILENAME: &'static str = "train/start";
const PFILENAME: &'static str = "train/stop";
const S1FILENAME: &'static str = "train/stop1";
const P1FILENAME: &'static str = "train/start1";
const DFILENAME: &'static str = "train/pwm";

const CODON: &'static [char] = &['A', 'C', 'G', 'T', 'N'];
const CODON_INDEL: &'static [char] = &['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', 'x'];

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

#[derive(Clone, Debug)]
pub struct HMM {
  pub initial_state: Vec<f64>,

  pub tr: Vec<f64>,

  pub e_m_1: Vec<Vec<Vec<f64>>>,
  pub e_m: Vec<Vec<Vec<f64>>>,

  pub tr_r_r: Vec<Vec<f64>>,
  pub tr_i_i: Vec<Vec<f64>>,
  pub tr_m_i: Vec<Vec<f64>>,

  pub tr_s: Vec<Vec<f64>>,
  pub tr_e: Vec<Vec<f64>>,
  pub tr_s_1: Vec<Vec<f64>>,
  pub tr_e_1: Vec<Vec<f64>>,

  pub s_dist: Vec<f64>,
  pub e_dist: Vec<f64>,
  pub s1_dist: Vec<f64>,
  pub e1_dist: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct Train {
  pub trans: Vec<Vec<Vec<Vec<f64>>>>,
  pub rtrans: Vec<Vec<Vec<Vec<f64>>>>,

  pub noncoding: Vec<Vec<Vec<f64>>>,
  pub start: Vec<Vec<Vec<f64>>>,
  pub stop: Vec<Vec<Vec<f64>>>,
  pub start1: Vec<Vec<Vec<f64>>>,
  pub stop1: Vec<Vec<Vec<f64>>>,

  pub s_dist: Vec<Vec<f64>>,
  pub e_dist: Vec<Vec<f64>>,
  pub s1_dist: Vec<Vec<f64>>,
  pub e1_dist: Vec<Vec<f64>>,
}
impl Train {
  // This is a static method
  // Static methods don't need to be called by an instance
  // These methods are generally used as constructors
  pub fn from_file() -> Train {
    let mut result = Train {
      trans: vec![vec![vec![vec![0.0; 4]; 16]; 6]; 44],
      rtrans: vec![vec![vec![vec![0.0; 4]; 16]; 6]; 44],

      noncoding: vec![vec![vec![0.0; 4]; 4]; 44],
      start: vec![vec![vec![0.0; 64]; 61]; 44],
      stop: vec![vec![vec![0.0; 64]; 61]; 44],
      start1: vec![vec![vec![0.0; 64]; 61]; 44],
      stop1: vec![vec![vec![0.0; 64]; 61]; 44],

      s_dist: vec![vec![0.0; 6]; 44],
      e_dist: vec![vec![0.0; 6]; 44],
      s1_dist: vec![vec![0.0; 6]; 44],
      e1_dist: vec![vec![0.0; 6]; 44],
    };
    result.load_m_state();
    result.load_m_1_state();
    result.load_noncoding_state();
    result.load_start_state();
    result.load_stop_state();
    result.load_start_1_state();
    result.load_stop_1_state();
    result.load_pwm_dist();
    result
  }
  fn load_m_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/gene.";
    let file = File::open(MFILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for i in 0..6 {
        for j in 0..16 {
          let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
          let tup: (f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
          self.trans[p][i][j][0] = tup.0.log2();
          self.trans[p][i][j][1] = tup.1.log2();
          self.trans[p][i][j][2] = tup.2.log2();
          self.trans[p][i][j][3] = tup.3.log2();
        }
      }
    }
  }

  fn load_m_1_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/rgene.";
    let file = File::open(M1FILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for i in 0..6 {
        for j in 0..16 {
          let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
          let tup: (f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
          self.rtrans[p][i][j][0] = tup.0.log2();
          self.rtrans[p][i][j][1] = tup.1.log2();
          self.rtrans[p][i][j][2] = tup.2.log2();
          self.rtrans[p][i][j][3] = tup.3.log2();
        }
      }
    }
  }

  fn load_noncoding_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/noncoding.";
    let file = File::open(NFILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..4 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        let tup: (f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
        self.noncoding[p][j][0] = tup.0.log2();
        self.noncoding[p][j][1] = tup.1.log2();
        self.noncoding[p][j][2] = tup.2.log2();
        self.noncoding[p][j][3] = tup.3.log2();
      }
    }
  }

  fn load_pwm_dist(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/noncoding.";
    let file = File::open(DFILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header

      let dists = [
        &mut self.s_dist,
        &mut self.e_dist,
        &mut self.s1_dist,
        &mut self.e1_dist,
      ];
      for i in 0..4 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        let tup: (f64, f64, f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
        dists[i][p][0] = tup.0;
        dists[i][p][1] = tup.1;
        dists[i][p][2] = tup.2;
        dists[i][p][3] = tup.3;
        dists[i][p][4] = tup.4;
        dists[i][p][5] = tup.5;
      }
    }
  }

  fn load_start_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/start.";
    let file = File::open(SFILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.start[p][j][k] = prob.log2();
        }
      }
    }
  }

  fn load_stop_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/stop.";
    let file = File::open(PFILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.stop[p][j][k] = prob.log2();
        }
      }
    }
  }

  fn load_start_1_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/stop1.";
    let file = File::open(S1FILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.start1[p][j][k] = prob.log2();
        }
      }
    }
  }

  fn load_stop_1_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/start1.";
    let file = File::open(P1FILENAME).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.stop1[p][j][k] = prob.log2();
        }
      }
    }
  }
}

impl HMM {
  pub fn from_file(train_file: &str) -> HMM {
    let mut result = HMM {
      initial_state: vec![0.0; 29],

      tr: vec![0.0; 14],

      e_m_1: vec![vec![vec![0.0; 4]; 16]; 6],
      e_m: vec![vec![vec![0.0; 4]; 16]; 6],

      tr_r_r: vec![vec![0.0; 4]; 4],
      tr_i_i: vec![vec![0.0; 4]; 4],
      tr_m_i: vec![vec![0.0; 4]; 4],

      tr_s: vec![vec![0.0; 64]; 61],
      tr_e: vec![vec![0.0; 64]; 61],
      tr_s_1: vec![vec![0.0; 64]; 61],
      tr_e_1: vec![vec![0.0; 64]; 61],

      s_dist: vec![0.0; 6],
      e_dist: vec![0.0; 6],
      s1_dist: vec![0.0; 6],
      e1_dist: vec![0.0; 6],
    };
    result.load_transition(&train_file);
    result
  }

  fn load_transition(&mut self, train_file: &str) {
    let file = File::open(train_file).expect("An error occured while opening the training file.");
    let reader = BufReader::new(file);

    for (index, line) in reader.lines().enumerate() {
      let line = line.expect("Something went wrong while reading the training file.");
      /* transition */
      if index > 0 && index < 15 {
        let tup: (String, f64) =
          parse_string(&line).expect("Unable to process 'Transition' from training file.");
        let tr = tr2int(&tup.0);
        self.tr[tr] = tup.1.log2();
      }
      /* transition MI */
      else if index > 15 && index < 32 {
        let tup: (char, char, f64) =
          parse_string(&line).expect("Unable to process 'TransitionMI' from training file.");
        let start = nt2int(tup.0);
        let end = nt2int(tup.1);
        self.tr_m_i[start][end] = tup.2.log2();
      }
      /* transition II */
      else if index > 32 && index < 49 {
        let tup: (char, char, f64) =
          parse_string(&line).expect("Unable to process 'TransitionII' from training file.");
        let start = nt2int(tup.0);
        let end = nt2int(tup.1);
        self.tr_i_i[start][end] = tup.2.log2();
      }
      /* PI */
      else if index > 49 {
        let tup: (String, f64) =
          parse_string(&line).expect("Unable to process 'PI' from training file.");
        self.initial_state[index - 50] = tup.1.log2();
      }
    }
  }
}

/*
Code for loading the relevant information into the HMM struct
 */

pub fn get_prob_from_cg(hmm: &mut HMM, train: &Train, seq: &String) -> usize {
  //change from void to int, Ye, April 18, 2016
  let mut cg_count: usize = 0;
  let len_seq = seq.len();

  for c in seq.chars() {
    cg_count += match c {
      'c' | 'C' | 'g' | 'G' => 1,
      _ => 0,
    }
  }
  let cg_count_i = ((((cg_count as f64 * 1.0) / len_seq as f64) * 100.0).floor() as i32) - 26;
  if cg_count_i < 0 {
    cg_count = 0;
  } else if cg_count_i > 43 {
    cg_count = 43;
  }
  else {
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

// TODO: move these helpers to separate file.
/*
HELPERS
*/

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

fn tr2int(tr: &String) -> usize {
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

fn nt2int_rc_indel(nt: &char) -> usize {
  match nt {
    'A' => 3,
    'C' => 2,
    'G' => 1,
    'T' => 0,
    'a' => 8,
    'c' => 7,
    'g' => 6,
    't' => 5,
    'n' => 9,
    'x' => 10,
    _ => 4,
  }
}

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

pub fn get_rc_dna(dna: &Vec<char>) -> Vec<char> {
  let len = dna.len();
  let mut result = vec!['\0'; len];
  for (i, c) in dna.iter().enumerate() {
    result[len - i - 1] = CODON[nt2int_rc(&*c)]
  }
  result
}

pub fn get_rc_dna_indel(dna: &Vec<char>) -> Vec<char> {
  let len = dna.len();
  let mut result = vec!['\0'; len];
  for (i, c) in dna.iter().enumerate() {
    result[len - i - 1] = CODON_INDEL[nt2int_rc_indel(&*c)]
  }
  result
}

fn strlen(s: &Vec<char>) -> usize { // TODO: cleanup
	let mut result = 0;
	for c in s.iter() {
		if *c != '\0' {
			result += 1;
    }
    else{
      break;
    }
	}
	result
}

pub fn get_protein(dna: &Vec<char>, strand: bool, wholegenome: bool) -> Vec<char> {
  let len = strlen(dna);
  let mut protein = vec![' '; len / 3];

  if strand {
    for i in (0..len).step_by(3) {
      if i/3 < len/3 {
      protein[i / 3] = CODON_CODE[trinucleotide_pep(&dna[i], &dna[i + 1], &dna[i + 2])];
      }
    }
  } else {
    for i in (0..len).step_by(3) {
      if (len - i) / 3 > 0 {
      protein[(len - i) / 3 - 1] =
        ANTI_CODON_CODE[trinucleotide_pep(&dna[i], &dna[i + 1], &dna[i + 2])];
      }
    }
  }

  if protein[len / 3 - 1] == '*' {
    //remove the ending *
    protein[len / 3 - 1] = '\0';
  }

  if wholegenome {
    return protein; //short reads, skip
  }

  if strand {
    let s = trinucleotide_pep(&dna[0], &dna[1], &dna[2]);
    if s == trinucleotide_pep(&'G', &'T', &'G') || s == trinucleotide_pep(&'T', &'T', &'G') {
      protein[0] = 'M';
    }
  } else {
    let s = trinucleotide_pep(&dna[len - 3], &dna[len - 2], &dna[len - 1]);
    if s == trinucleotide_pep(&'C', &'A', &'C') || s == trinucleotide_pep(&'C', &'A', &'A') {
      protein[0] = 'M';
    }
  }
  protein
}

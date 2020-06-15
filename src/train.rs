use super::dna_helpers::{nt2int, tr2int};
use super::helpers::get_executable_path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use whiteread::parse_string;

/**
 * train.rs
 * ========
 * This file contains code responsible for reading training files and building the respective Train and HMM structs
 * from those files.
 */

const MFILENAME: &'static str = "train/gene";
const M1FILENAME: &'static str = "train/rgene";
const NFILENAME: &'static str = "train/noncoding";
const SFILENAME: &'static str = "train/start";
const PFILENAME: &'static str = "train/stop";
const S1FILENAME: &'static str = "train/stop1";
const P1FILENAME: &'static str = "train/start1";
const DFILENAME: &'static str = "train/pwm";

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
  /**
   * This method will build a new Train struct from the train/ folder stored in the same path as the
   * executable.
   */
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
  /**
   * This method wil load the train/gene file into the trans field of our struct
   */
  fn load_m_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/gene.";
    let mut path = get_executable_path();
    path.push(MFILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for i in 0..6 {
        for j in 0..16 {
          let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
          let tup: (f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
          self.trans[p][i][j][0] = tup.0.ln();
          self.trans[p][i][j][1] = tup.1.ln();
          self.trans[p][i][j][2] = tup.2.ln();
          self.trans[p][i][j][3] = tup.3.ln();
        }
      }
    }
  }
  /**
   * This method wil load the train/rgene file into the rtrans field of our struct
   */
  fn load_m_1_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/rgene.";
    let mut path = get_executable_path();
    path.push(M1FILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for i in 0..6 {
        for j in 0..16 {
          let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
          let tup: (f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
          self.rtrans[p][i][j][0] = tup.0.ln();
          self.rtrans[p][i][j][1] = tup.1.ln();
          self.rtrans[p][i][j][2] = tup.2.ln();
          self.rtrans[p][i][j][3] = tup.3.ln();
        }
      }
    }
  }

  /**
   * This method wil load the train/noncoding file into the noncoding field of our struct
   */
  fn load_noncoding_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/noncoding.";
    let mut path = get_executable_path();
    path.push(NFILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..4 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        let tup: (f64, f64, f64, f64) = parse_string(&line).expect(READ_ERROR);
        self.noncoding[p][j][0] = tup.0.ln();
        self.noncoding[p][j][1] = tup.1.ln();
        self.noncoding[p][j][2] = tup.2.ln();
        self.noncoding[p][j][3] = tup.3.ln();
      }
    }
  }

  /**
   * This method wil load the train/pwm file into the s_dist/e_dist/s1_dist/e1_dist fields of our struct
   */
  fn load_pwm_dist(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/pwm.";
    let mut path = get_executable_path();
    path.push(DFILENAME);
    let file = File::open(path).expect(READ_ERROR);
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
  /**
   * This method wil load the train/start file into the start field of our struct
   */
  fn load_start_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/start.";
    let mut path = get_executable_path();
    path.push(SFILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.start[p][j][k] = prob.ln();
        }
      }
    }
  }

  /**
   * This method wil load the train/stop file into the stop field of our struct
   */
  fn load_stop_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/stop.";
    let mut path = get_executable_path();
    path.push(PFILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.stop[p][j][k] = prob.ln();
        }
      }
    }
  }

  /**
   * This method wil load the train/start file into the start1 fields of our struct
   */
  fn load_start_1_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/stop1.";
    let mut path = get_executable_path();
    path.push(S1FILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.start1[p][j][k] = prob.ln();
        }
      }
    }
  }
  /**
   * This method wil load the train/stop file into the stop1 fields of our struct
   */
  fn load_stop_1_state(&mut self) {
    const READ_ERROR: &'static str = "Something went wrong while reading train/start1.";
    let mut path = get_executable_path();
    path.push(P1FILENAME);
    let file = File::open(path).expect(READ_ERROR);
    let mut lines = BufReader::new(file).lines();

    for p in 0..44 {
      lines.next(); // first line is header
      for j in 0..61 {
        let line = lines.next().expect(READ_ERROR).expect(READ_ERROR);
        for (k, value) in line.split_whitespace().enumerate() {
          assert!(k < 64, READ_ERROR);
          let prob = value.parse::<f64>().expect(READ_ERROR);
          self.stop1[p][j][k] = prob.ln();
        }
      }
    }
  }
}

impl HMM {
  /**
   * This method will yield the HMM struct based on the specified train file path
   */
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

  /**
   * This method will read the contents of the specified training file path into the HMM struct
   */
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
        self.tr[tr] = tup.1.ln();
      }
      /* transition MI */
      else if index > 15 && index < 32 {
        let tup: (char, char, f64) =
          parse_string(&line).expect("Unable to process 'TransitionMI' from training file.");
        let start = nt2int(tup.0);
        let end = nt2int(tup.1);
        self.tr_m_i[start][end] = tup.2.ln();
      }
      /* transition II */
      else if index > 32 && index < 49 {
        let tup: (char, char, f64) =
          parse_string(&line).expect("Unable to process 'TransitionII' from training file.");
        let start = nt2int(tup.0);
        let end = nt2int(tup.1);
        self.tr_i_i[start][end] = tup.2.ln();
      }
      /* PI */
      else if index > 49 {
        let tup: (String, f64) =
          parse_string(&line).expect("Unable to process 'PI' from training file.");
        self.initial_state[index - 50] = tup.1.ln();
      }
    }
  }
}

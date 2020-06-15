use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

/**
 * helpers.rs
 * ==============
 * This file contains general helper methods related to file I/O,
 * and String manipulations.
 */

/**
 * This method will swap the error message for file opening/creation
 */
pub fn create_file_if_not_exists(path: &str) -> File {
  match Path::new(path).exists() {
    true => File::create(path).expect("Error: unable to open output file"),
    false => File::create(path).expect("Error: unable to create output file"),
  }
}

/*
* Helper method to catch any possible errors in writing to file.
*/
pub fn write_data(output: &mut File, data: String) {
  match write!(output, "{}", data) {
    Err(e) => eprintln!("Error: {:?}", e,),
    _ => (),
  }
}

/*
* Helper method to get executable path (e.g. used to find relative location of the train/ foder)
*/
pub fn get_executable_path() -> PathBuf {
  let path = match env::current_exe() {
    Ok(mut path) => {
      path.pop();
      path
    }
    Err(_) => panic!("[Error] Current executable path does not exist"),
  };
  path
}

/**
 * Helper method to emulate the behaviour of strlen
 */
pub fn strlen(s: &String) -> usize {
  let mut result = 0;
  for c in s.chars() {
    if c != '\0' {
      result += 1;
    } else {
      break;
    }
  }
  result
}
/**
 * Helper method to emulate the behaviour of strncpy
 */
pub fn strncpy(target: &mut Vec<char>, source: &Vec<char>, start_pos: usize, len: usize) {
  let sourcelen = source.len();
  for i in start_pos..(start_pos + len) {
    if i < sourcelen {
      target[i - start_pos] = source[i];
    } else {
      target[i - start_pos] = '\0';
    }
  }
}

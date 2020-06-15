use std::fs::File;
use std::path::Path;
use std::io::{Write};
use std::path::PathBuf;
use std::env;



pub fn create_file_if_not_exists(path: &str) -> File {
    match Path::new(path).exists() {
        true => File::create(path).expect("Error: unable to open output file"),
        false => File::create(path).expect("Error: unable to create output file"),
    }
}

pub fn write_data(output: &mut File, data: String) {
    match write!(output, "{}", data) {
        Err(e) => eprintln!("Error: {:?}", e,),
        _ => (),
    }
}

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
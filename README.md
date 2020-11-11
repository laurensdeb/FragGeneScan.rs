# FragGeneScan.rs: Rust Implementation of FragGeneScan 
This repository contains a Rust implementation of the popular FragGeneScan utility. It fixes many common race conditions and parallellism issues with the original FragGeneScan-tool and improves upon code quality.

## Technologies used

* Rust
* Rayon
* Rust-Bio
* Clap
* Whiteread (for parsing the FragGeneScan input files)

## Usage
To build the project:
```sh
cargo build
```

To make a release build of the project:
```sh
cargo build --release
```

Get the utility's help documentation
```sh
./fgsrs --help
````

Process one of the example files (output of FragGeneScan 1.31 is also present in the `/example` folder):
```sh
./fgsrs -w 0 -p 20 -t train/454_10 < example/NC_000913-454.fna > NC_000913-454-fgs.faa 
```

## License
FragGeneScan.rs is made available under the MIT License. Reference to the author and this repository however is appreciated.

## References
Based upon the original implementation:

`Mina Rho, Haixu Tang, and Yuzhen Ye. FragGeneScan: Predicting Genes in Short and Error-prone Reads. Nucl. Acids Res., 2010 doi: 10.1093/nar/gkq747`
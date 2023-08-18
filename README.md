# GenomeTK

This project that was created in my free time to get familiarized with developing
command line interfaces in Rust and is not supposed to be used for serious projects.
The Rust programming concepts that were covered in this project were:
    1. Writing Rust code
    2. Using external crates
    3. Creating multi-module applications
    4. Creating command-line arguments
    5. Testing

## Credit where credit is due

This project was inspired by the [seqtk](https://github.com/lh3/seqtk)
project which is a complete toolkit for handling biological sequences through
the command line.

## Specifications

This application consists of three modules:

 1. Generate
 The generate module generates a specified number of random sequences with
 specified length. If mutation rate is provided the application will
 create a seed sequence and mutate it

 2. Tabulate
 This module accepts a fasta file as input and create a tab-delimited file
 with basic information about the sequences such an a example table is shown below:
    |SeqName|Length|GC%|A|T|G|C|N|A%|T%|G%|C%|N%|
    |-------|------|---|-|-|-|-|-|--|--|--|--|--|
    |Seq1|100|50|25|25|25|25|0|25|25|25|25|0|
    |...|.|.|.|.|.|.|.|.|.|.|.|.|.|.|
    |...|.|.|.|.|.|.|.|.|.|.|.|.|.|.|
    |SeqN|100|40|30|30|20|20|0|30|30|20|20|0|

 3. Subsample
 This module subsamples a specified number of sequences from a fasta file.
 By providing a file with a list of sequence identifiers (one per line) the
 module extracts these specific sequences to a new file.

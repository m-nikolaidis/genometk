extern crate bio;
use bio::io::fasta;
use std::fs::File;
use std::io::{ BufWriter, BufReader, BufRead};

fn calculate_fasta_records(inputfile: &str) -> i64 {
    let mut records = fasta::Reader::from_file(&inputfile).unwrap().records();
    let mut counter = 0;
    while let Some(Ok(_)) = records.next() {
        counter += 1;
    }
    return counter;
}

pub fn subsample_sequences(input: &str, output: &str, num_sequences: usize){
    // Random sample of size num_sequences from input fasta file
    fastrand::seed(0);
    let total_sequences: i64 = calculate_fasta_records(&input);
    if num_sequences > total_sequences as usize {
        panic!("The number of sequences requested for the sample is greater than the total number of sequences in the input file");
    }
    let record_ids: Vec<i64> = fastrand::choose_multiple(1..total_sequences, num_sequences);

    // Open output file
    let outfile_res = File::create(output);
    let outfile = match outfile_res {
        Ok(f) => f,
        Err(_) => panic!("Error creating output file"),
    };

    let mut writer = fasta::Writer::new(BufWriter::new(outfile));
    let mut records = fasta::Reader::from_file(&input).unwrap().records();
    let mut counter = 0;
    while let Some(Ok(record)) = records.next() {
        if record_ids.contains(&counter) {
            writer.write_record(&record).unwrap();
        }
        counter += 1;
    }
}

fn read_lines(f: File) -> Vec<String> {
        let list_reader = BufReader::new(f);
        let list: Vec<String> = list_reader.lines().map(|line| {
            line.unwrap().trim().to_string()
        }).collect();
    return list;
}

pub fn subsample_list(input: &str, output: &str, list_file: &str){
    let file_reader = File::open(list_file);
    if let Err(_) = file_reader {
        panic!("Could not open {}", input);
    }
    let seqlist = read_lines(file_reader.unwrap());

    // Open output file
    let outfile_res = File::create(output);
    let outfile = match outfile_res {
        Ok(f) => f,
        Err(_) => panic!("Error creating output file"),
    };

    // Read fasta file and stream to output
    let mut writer = fasta::Writer::new(BufWriter::new(outfile));
    let mut records = fasta::Reader::from_file(&input).unwrap().records();
    while let Some(Ok(record)) = records.next(){
        if seqlist.contains(&record.id().to_string()) {
            writer.write_record(&record).unwrap();
        }
    }
}

#[test]
fn test_subsample_sequences() {
    let input = "test_files/test.fasta";
    let output = "test_files/test_subsample.fasta";
    let num_sequences = 2;
    subsample_sequences(&input, &output, num_sequences);
    let mut records = fasta::Reader::from_file(&output).unwrap().records();
    let mut counter = 0;
    while let Some(Ok(_)) = records.next() {
        counter += 1;
    }
    assert_eq!(counter, num_sequences);
}

#[test]
fn test_subsample_sequences_more_than_max(){
    let input = "test_files/test.fasta";
    let output = "test_files/test_subsample.fasta";
    let num_sequences = 100;
    let result = std::panic::catch_unwind(|| subsample_sequences(&input, &output, num_sequences));
    assert!(result.is_err());
}

#[test]
fn test_subsample_list() {
    let input = "test_files/test.fasta";
    let output = "test_files/test_subsample_list.fasta";
    let list_file = "test_files/test_list.txt";
    subsample_list(&input, &output, &list_file);
    let records = fasta::Reader::from_file(&output).unwrap().records();
    let output_records: Vec<String> = records.map(|rec| rec.unwrap().id().to_string()).collect();
    let desired_records: Vec<String> = vec!["Seq1".to_string(), "Seq3".to_string()];
    assert_eq!(output_records, desired_records);
}

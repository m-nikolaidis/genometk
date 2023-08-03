extern crate bio;
use bio::io::fasta;
use std::fs::File;
use std::io::{BufWriter, Write};

fn gather_length(record: &fasta::Record) -> usize {
    record.seq().len()
}

fn calculate_gc(g: usize, c: usize) -> usize {
    return g + c;
}

fn gather_freq(record: &fasta::Record) -> [usize; 5] {
    let mut a = 0;
    let mut t = 0;
    let mut g = 0;
    let mut c = 0;
    let mut n = 0;
    for nucl in record.seq().iter() {
        if *nucl == b'A' {
            a += 1;
        }
        if *nucl == b'T' {
            t += 1;
        }
        if *nucl == b'G' {
            g += 1;
        }
        if *nucl == b'C' {
            c += 1;
        }
        if *nucl == b'N' {
            n += 1;
        }
    }
    return [a, t, g, c, n];
}

fn gather_info(record: &fasta::Record) -> (usize, usize, [usize; 5]) {
    let length = gather_length(&record);
    let freq = gather_freq(&record);
    let gc = calculate_gc(freq[2], freq[3]);
    return (length, gc, freq);
}

pub fn tabulate(inputfile: &str, out_filename: &str) {

    // Create output file
    let outfile_res = File::create(out_filename);
    let outfile = match outfile_res {
        Ok(f) => f,
        Err(_) => panic!("Error creating output file"),
    };
    let mut writer = BufWriter::new(outfile);

    // Header
    writeln!(writer, "ID\tLength\tGC\tA\tT\tG\tC\tN").unwrap();

    // Read input fasta file and write
    let mut records = fasta::Reader::from_file(&inputfile).unwrap().records();
    // This unwrap does not gracefully handle the error

    while let Some(Ok(record)) = records.next() {
        let info = gather_info(&record);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.id(),
            info.0,
            info.1,
            info.2[0],
            info.2[1],
            info.2[2],
            info.2[3],
            info.2[4],
        )
        .unwrap();
    }
}

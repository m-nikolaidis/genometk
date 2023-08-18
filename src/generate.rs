use bio::io::fasta;
use rayon::prelude::*;
use std::fs::{File, remove_file};
use std::io::BufWriter;
use std::sync::{Arc, Mutex};
use fastrand;


fn create_seq(seqlen: i64, fastrand_base: &mut fastrand::Rng) -> Vec<u8> {
    let nucleotides = [b'A', b'C', b'T', b'G'];
    (0..seqlen)
        .map(|_| nucleotides[fastrand_base.usize(0..4)])
        .collect::<Vec<u8>>()
}

fn mutate_reference(refseq: &Vec<u8>, mut_rate: &f64, nucleotides: &[u8; 4]) -> Vec<u8> {
    let seq = refseq
        .clone()
        .into_iter()
        .map(|nucl| {
            if fastrand::f64() < *mut_rate {
                return nucleotides[fastrand::usize(0..4)];
            }
            return nucl;
        })
        .collect::<Vec<u8>>();
    return seq;
}

pub fn generate_sequences(outfname: &str, seqlen: i64, maxseqs: i64, mut_rate: Option<&f64>) {

    // Set up fastrand generator
    let mut fastrand_base = fastrand::Rng::new();
    fastrand_base.seed(1234);

    // Set up nucleotides
    let nucleotides: [u8; 4] = [b'A', b'C', b'T', b'G'];

    // Set up buffered writer
    let outfile_res = File::create(outfname);
    let outfile = match outfile_res {
        Ok(f) => f,
        Err(_) => panic!("Error creating output file"),
    };
    let handle = BufWriter::new(outfile);
    let writer = Arc::new(Mutex::new(fasta::Writer::from_bufwriter(handle)));

    let refseq = create_seq(seqlen, &mut fastrand_base);
    if mut_rate.is_some() {
        let mut_rate = mut_rate.unwrap();
        let record = fasta::Record::with_attrs("Refseq", None, &refseq);
        writer.lock().unwrap().write_record(&record).unwrap();
        let _r = (0..maxseqs - 1)
            .into_par_iter()
            .map(|i| {
                let mutseq = mutate_reference(&refseq, &mut_rate, &nucleotides);
                let seqname = format!("seq{}", i);
                let record = fasta::Record::with_attrs(&seqname, None, &mutseq);
                writer.lock().unwrap().write_record(&record).unwrap();
            })
            .collect::<Vec<_>>();
    } else {
        let _r = (0..maxseqs)
            .into_iter()
            .map(|i| {
                let mutseq = create_seq(seqlen, &mut fastrand_base);
                let seqname = format!("seq{}", i);
                let record = fasta::Record::with_attrs(&seqname, None, &mutseq);
                writer.lock().unwrap().write_record(&record).unwrap();
            })
            .collect::<Vec<_>>();
    }
}

#[test]
fn test_generate_for_number_of_sequences() {
    let outfname = "test.fasta";
    let seqlen = 10;
    let maxseqs = 5;
    generate_sequences(outfname, seqlen, maxseqs, None);
    let reader = fasta::Reader::from_file(outfname).unwrap();
    let mut count = 0;
    let mut len = 0;
    for record in reader.records() {
        let record = record.unwrap();
        count += 1;
        let tmp_len = record.seq().len();
    
        if tmp_len > len{
            len = tmp_len;
        }
    }
    let len = len as i64;
    assert_eq!(count, maxseqs);
    assert_eq!(len, seqlen);
    remove_file(outfname).unwrap();
}

#[test]
#[cfg(debug_assertions)]
fn test_generate_no_mut_rate() {
    let outfname = "test2.fasta";
    let seqlen = 10;
    let maxseqs = 100;
    generate_sequences(outfname, seqlen, maxseqs, None);
    let reader = fasta::Reader::from_file(outfname).unwrap();
    let mut count = 0;
    let records = reader.records().map(|r| r.ok().unwrap()).collect::<Vec<_>>();
    let ref_seq = &records[0].seq();
    let mut total_changes = 0;
    for record in &records {
        count += 1;
        if count == 1 {
            continue;
        }
        let mut record_changes = 0;
        for nucl_pos in 0..record.seq().len() {
            if record.seq()[nucl_pos] != ref_seq[nucl_pos]{
            record_changes += 1;
            }
        }
        if record_changes > total_changes{
            total_changes = record_changes;
        }
    }
    debug_assert!(total_changes > 1);
    remove_file(outfname).unwrap();
}

#[test]
#[cfg(debug_assertions)]
fn test_generate_mut_rate() {
    let outfname = "test3.fasta";
    let seqlen = 10_000;
    let maxseqs = 100;
    let mut_rate = 1e-6;
    generate_sequences(outfname, seqlen, maxseqs, Some(&mut_rate));
    let reader = fasta::Reader::from_file(outfname).unwrap();
    let mut count = 0;
    let records = reader.records().map(|r| r.ok().unwrap()).collect::<Vec<_>>();
    let ref_seq = &records[0].seq();
    let mut total_changes = 0;
    for record in &records {
        count += 1;
        if count == 1 {
            continue;
        }
        let mut record_changes = 0;
        for nucl_pos in 0..record.seq().len() {
            if record.seq()[nucl_pos] != ref_seq[nucl_pos]{
            record_changes += 1;
            }
        }
        if record_changes > total_changes{
            total_changes = record_changes;
        }
    }
    debug_assert!(total_changes <= 1);
    remove_file(outfname).unwrap();
}

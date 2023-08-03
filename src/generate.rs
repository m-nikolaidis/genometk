use bio::io::fasta;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::sync::{Arc, Mutex};


fn create_seq(seqlen: i64) -> Vec<u8> {
    let nucleotides = [b'A', b'C', b'T', b'G'];
    (0..seqlen)
        .map(|_| nucleotides[fastrand::usize(0..4)])
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

    let refseq = create_seq(seqlen);
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
            .into_par_iter()
            .map(|i| {
                let mutseq = create_seq(seqlen);
                let seqname = format!("seq{}", i);
                let record = fasta::Record::with_attrs(&seqname, None, &mutseq);
                writer.lock().unwrap().write_record(&record).unwrap();
            })
            .collect::<Vec<_>>();
    }
}

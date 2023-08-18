// use camino::Utf8PathBuf;
use clap::{arg, Command, ArgGroup};
mod generate;
mod tabulate;
mod subsample;

fn cli() -> Command {
    let cmd = Command::new("genometk")
        .bin_name("genometk")
        .version("0.1.0")
        .author("M. Nikolaidis")
        .about("A tool for generating and manipulating genomes")
        .subcommand(
            Command::new("generate")
            .about("Generate a genome")
            .arg(arg!(-o <output> "The Output file").required(true))
            .arg(arg!(-l <seqlen> "The length of the sequences").required(true).value_parser(clap::value_parser!(i64)))
            .arg(arg!(-n <num_sequences> "The number of the final sequences").required(true).value_parser(clap::value_parser!(i64)))
            .arg(arg!(-m <mut_rate> "Mutation rate").required(false).value_parser(clap::value_parser!(f64)))
        )
        .subcommand(
            Command::new("tabulate")
                .about("Create table with information from genomes")
                .arg(arg!(-i <input> "Input fasta file").required(true))
                .arg(arg!(-o <output> "The Output file").required(true))
        )
        .subcommand(
            Command::new("sample")
                .about("Subsample a fasta file")
                .arg(arg!(-i <input> "Input fasta file").required(true))
                .arg(arg!(-o <output> "The Output file").required(true))
                .arg(arg!(-n <num_sequences> "The number of the final sequences").required(false).value_parser(clap::value_parser!(usize)))
                .arg(arg!(-f <list_file> "A file with the list of sequences to keep").required(false).value_parser(clap::value_parser!(String)))
                .group(ArgGroup::new("sample_type")
                    .args(&["num_sequences", "list_file"])
                    .required(true)
                )
        );
    return cmd;
}


fn main(){
    let matches = cli().get_matches();

    match matches.subcommand() {
        Some(("generate", subcmd_gen)) => {
            let seq_len = subcmd_gen.get_one::<i64>("seqlen").unwrap();
            let output = subcmd_gen.get_one::<String>("output").unwrap();
            let num_sequences = subcmd_gen.get_one::<i64>("num_sequences").unwrap();
            let mut_rate = subcmd_gen.get_one::<f64>("mut_rate");
            generate::generate_sequences(&*output, *seq_len, *num_sequences, mut_rate);
        },
        Some(("tabulate", subcmd_tab)) => {
            let input = subcmd_tab.get_one::<String>("input").unwrap().trim();
            let output = subcmd_tab.get_one::<String>("output").unwrap().trim();
            tabulate::tabulate(&input, &output);
        },
        Some(("sample", subcmd_sample)) => {
            let input = subcmd_sample.get_one::<String>("input").unwrap().trim();
            let output = subcmd_sample.get_one::<String>("output").unwrap().trim();
            let num_sequences = subcmd_sample.get_one::<usize>("num_sequences");
            let list_file = subcmd_sample.get_one::<String>("list_file");
            match num_sequences {
                Some(num) => subsample::subsample_sequences(&input, &output, *num),
                None => (),
            }
            match list_file {
                Some(lf) => subsample::subsample_list(&input, &output, &lf.trim()),
                None => (),
            }
        }
        _ => println!("No subcommand was used. Use the -h flag to see the available subcommands"),
    }

}

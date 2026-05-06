use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <type: fastq|vcf> <output_file>", args[0]);
        std::process::exit(1);
    }

    let mode = &args[1];
    let output = &args[2];
    
    let file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, file);

    if mode == "fastq" {
        // Generate 1M reads
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        for i in 0..1_000_000 {
            writeln!(writer, "@synthetic_read_{}\n{}\n+\n{}", i, seq, qual)?;
        }
    } else if mode == "vcf" {
        // Generate 100K variants
        writeln!(writer, "##fileformat=VCFv4.2")?;
        writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")?;
        for i in 1..=100_000 {
            let alt = if i % 10 == 0 { "A,C" } else if i % 5 == 0 { "ATG" } else { "T" };
            writeln!(writer, "chr1\t{}\t.\tA\t{}\t50\tPASS\tDP=100;AF=0.5\tGT:DP\t0/1:100", i * 10, alt)?;
        }
    } else {
        eprintln!("Unknown mode: {}", mode);
        std::process::exit(1);
    }

    Ok(())
}

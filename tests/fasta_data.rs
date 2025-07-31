// ntcomp: Sequencing data compression using SBWT and k-bounded matching statistics.
//
// Copyright 2025 Tommi Mäklin [tommi@maklin.fi].
//
// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.
//
// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
use random::Source;

fn random_nucleotide(
    rng: &mut random::Default
) -> u8 {
    match rng.read_u64() % 4 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => panic!("Impossible math")
    }
}

#[test]
fn random_fasta_data() {
    use std::io::Cursor;
    use std::io::Read;
    use std::io::Seek;
    use std::io::Write;
    use ntcomp::encode_sequence;
    use ntcomp::encode_file_header;
    use ntcomp::decode_file_header;

    let mut rng = random::default(20250731);
    let mut contigs: Vec<Vec<u8>> = Vec::new();

    let n_contigs = 7;
    let max_contig_len = 2000;

    for _ in 0..n_contigs {
        let contig_len = rng.read_u64() % max_contig_len;
        let contig: Vec<u8> = (0..contig_len).map(|_| random_nucleotide(&mut rng)).collect();
        contigs.push(contig);
    }

    let mut sbwt_build_options = kbo::BuildOpts::default();
    sbwt_build_options.k = 255;
    sbwt_build_options.dedup_batches = false;
    sbwt_build_options.temp_dir = None;
    sbwt_build_options.add_revcomp = true;
    sbwt_build_options.build_select = true;

    let (sbwt, lcs) = kbo::build(&contigs, sbwt_build_options);

    let n_bases = contigs.iter().map(|x| x.len()).sum::<usize>();
    let mut buf = Cursor::new(vec![0; n_bases + n_contigs*11]);
    let block_size = 3;

    let header_bytes = encode_file_header(0,0,0,0).unwrap();
    let _ = buf.write_all(&header_bytes);

    let mut u64_encoding: Vec<u64> = Vec::new();
    let mut num_records = 0;
    contigs.iter().for_each(|nucleotides| {
        num_records += 1;
        u64_encoding.append(&mut encode_sequence(&nucleotides, &sbwt, &lcs));

        if num_records % block_size == 0 {
            let _ = ntcomp::write_block_to(&u64_encoding, num_records, &mut buf);
            u64_encoding.clear();
        }
    });
    if !u64_encoding.is_empty() {
        let _ = ntcomp::write_block_to(&u64_encoding, num_records, &mut buf);
    }
    let _ = buf.rewind();

    // File header
    let mut header_bytes: [u8; 32] = [0_u8; 32];
    let _ = buf.read_exact(&mut header_bytes);
    let _file_header = decode_file_header(&header_bytes);

    let mut i = 0;
    let mut total_bases = 0;
    while let Ok(records) = ntcomp::decode_block(&_file_header, &sbwt, &mut buf) {
        records.iter().for_each(|nucleotides| {
            eprintln!("{},{},{}", nucleotides.len(), contigs[i].len(), i);
            assert_eq!(*nucleotides, contigs[i]);
            total_bases += nucleotides.len();
            i += 1;
        });
    }

    assert_eq!(n_contigs, i);
    assert_eq!(n_bases, total_bases);
}

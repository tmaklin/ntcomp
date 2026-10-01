# ntcomp
(Proof of concept) Reference-based microbial sequencing data compression using [SBWT](https://docs.rs/sbwt) and [_k_-bounded matching statistics](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1).

Implements two algorithms for compressing/decompressing microbial genome assemblies:
- nlz: Nucleotide relative Lempel-Ziv, encodes the sequence contents of an assembly using its SBWT as the dictionary.
- nlz-og: Nucleotide relative Lempel-Ziv on graphs, stores the nlz encodings for closely related assemblies in a graph that deduplicates shared segments.

*Note* this project is currently a work-in-progress and should only be used for research purposes.

## Installation
``` sh
git clone https://github.com/tmaklin/ntcomp
cd ntcomp
cargo build --release
```
The built binary is located at `target/release/ntcomp`.

## Usage
### Setup example data
``` sh
wget https://a3s.fi/maklinto-2006181-pub/ntcomp-test.tar
tar -xf ntcomp-test.tar
```

### Preparing an index
Build from a reference sequence
``` sh
ntcomp build -o index -k91 test/GCA_964037205.1_30348_1_60_genomic.fna.gz
```

Larger values of `-k` produce better compression ratios at the cost of larger .lcs file size.

### nlz
#### Encoding fastX data
Query fastX data against an index and write the encoding
``` sh
ntcomp encode --index index test/ERR10498075.fastq.gz > encoded.dat
```

Encoding data requires both the .sbwt and .lcs files.

This removes quality scores.

#### Decoding encoded fastX data
Retrieve the encoded sequences from an index

``` sh
ntcomp decode --index index encoded.dat > decoded.fasta
```

In theory decoding requires only the .sbwt file but the tool will not run without the .lcs file present.

Decoding will also work with a rebuilt SBWT, as the construction algorithm is deterministic.

#### Verify
Requires installing [seqtk](https://github.com/lh3/seqtk)
``` sh
seqtk seq -A test/ERR10498075.fastq.gz | sed 's/ERR[0-9]*[.]\([0-9]*\).*$/seq.\1/g' > expected.fasta
diff -s decoded.fasta expected.fasta
```

### nlz-og
#### Encoding a collection of fastX files
```sh
ntcomp graph -i index query_1.fna query_2.fna query_3.fna > graph_3.dat
```

#### Decoding an encoded collection

``` sh
ntcomp graph -di index --fasta graph_3.dat > queries_1-3.fasta
```

## About
Works by finding the longest common suffix matches between the _k_-mers in an input nucleotide sequence and an SBWT index and encoding the input as pairs of (`longest common suffix length` and `suffix sequence location in the SBWT`).

If the index is a close match to the sequencing reads, such as an assembly from the same reads, and the reads are reasonably accurate, many of the _k_-mers in the reads are redundant and looking them up from the index uses less storage compared to storing their entire sequence.

## License
ntcomp is dual-licensed under the [MIT](LICENSE-MIT) and [Apache 2.0](LICENSE-APACHE) licenses.


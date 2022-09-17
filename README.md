# EDS
Automata based Degenerate String and Elastic Degenerate String parser


## Documentation

```
cargo doc --open
```

## EDS format

Degenerate letters are between brackets.
Variants are separated by commas.

```
ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATATTGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGTTG{A,G}CTC{T,G}TGAGTGAGCTTGCGAGATA
```

| Symbol | State             |
|--------|-------------------|
| D      | Degenerate Letter |
| S      | Solid String      |

| Symbol | Char                                                         |
|--------|--------------------------------------------------------------|
| N      | Nucleotide                                                   |
| C      | Comma                                                        |
| O      | Open bracket / degenerate letter start / solid string end    |
| E      | Closing bracket / degenerate letter end / solid string start |
|        |                                                              |

## Tooling

- Generate Elastic Degenerate Text from fasta file and VCF https://github.com/urbanslug/aedso
- Generate random Degenerate Text https://github.com/urbanslug/simed
- Generate random Elastic Degenerate Text https://github.com/webmasterar/EDSRand

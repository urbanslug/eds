# EDS
EDS parser


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

- Generate EDS from fasta file and VCF https://github.com/urbanslug/aedso
- Generate random EDS https://github.com/urbanslug/simed

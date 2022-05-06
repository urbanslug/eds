# EDS

EDS parser

```
ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATATTGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGTTG{A,G}CTC{T,G}TGAGTGAGCTTGCGAGATA
```
Parse an EDS into a DAG.


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

![EDS-automata](./Figures/eds-automata.png)

Degenerate letters are between brackets.
Individual seeds are separated by commas.

Generate EDS from fasta file and VCF
https://github.com/urbanslug/aedso

Generate random EDS
https://github.com/webmasterar/EDSRand

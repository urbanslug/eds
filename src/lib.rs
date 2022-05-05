/*!
# EDS

EDT parser

Degenerate letters are between brackets.
Individual seeds are separated by commas.
Empty seeds are allowed
`AT{TCC,AG,C,}AA`

Generate EDS from fasta file and VCF
[https://github.com/urbanslug/aedso]()

Generate random EDS
[https://github.com/webmasterar/EDSRand]()

```
use eds::EDT;

let ed_string = "AT{TCC,C,}AA";
let edt = EDT::from_str(ed_string);

// ---------------
// Size and Length
// ---------------
assert_eq!(edt.size() as usize, 8);
assert_eq!(edt.length() as usize, 5);

// -----
// Edges
// -----
// TODO: add

// --------
// Indexing
// --------

// 0123456
// ATTCCAA
//   C**
//   ***


assert_eq!(edt[0], vec![b'A']);
assert_eq!(edt[2], vec![b'T', b'C',  b'*']);
assert_eq!(edt[edt.p() - 1], vec![b'A']);
```

 */

use std::collections::HashSet;
use std::fmt::{self, Display};
use std::fs;
use std::ops::Index;

// ---------
// Constants
// ----------
// Minimum number of chars it can hold without re-allocating
// assume small viral pan-genome (50k)
const EXPECTED_COLS: usize = 50_000; // expect no seq to be over 3000 bases
const EXPECTED_ROWS: usize = 5; // expect a max of 5 variations
const WILDCARD: u8 = b'*';
const PAD: bool = true;

/// A `[column, row]` index
pub type Coordinate = [usize; 2];

/// A letter in the EDT. It can either be degenerate or solid
#[derive(Debug, PartialEq, Copy, Clone)]
enum Letter {
    Solid,
    Degenerate,
}

/// A single ASCII character expected while parsing the EDT.
/// Expected characters are
/// one of four characters A, T, C or G for nucleotides
/// , a comma to separate degenerate letters
/// { an open curly bracket to start a degenerate letter
/// } a closed curly bracket to end a degenerate letter
#[derive(Debug, PartialEq, Copy, Clone)]
enum Char {
    Nucleotide,
    Comma,
    Open,
    Close,
}

/// 0 => from (incoming nodes)
/// 1 => to (outgoing nodes)
#[derive(Debug)]
pub struct Edges(HashSet<Coordinate>, HashSet<Coordinate>);

impl Edges {
    #[allow(dead_code)]
    fn new() -> Self {
        Edges(HashSet::<Coordinate>::new(), HashSet::<Coordinate>::new())
    }

    // outgoing edges
    #[allow(dead_code)]
    fn from(&self) -> &HashSet<Coordinate> {
        &self.1
    }

    // incoming edges
    #[allow(dead_code)]
    fn to(&self) -> &HashSet<Coordinate> {
        &self.0
    }
}

/// The underlying Elastic Degenerate Text
pub struct EDT {
    pub data: Vec<Vec<u8>>,

    pub edges: Edges,

    /// deepest z
    z: usize,

    /// diameter of the graph or length of longest string in the possibility set
    p: usize,

    /// size (N)
    length: usize,

    /// length (n)
    size: usize,

    /// An array of sets of start indices
    /// Start indices for each nucleotide
    /// Index of where each char starts
    /// A = 0 T = 1 C = 2 G = 3
    start_indices: [HashSet<usize>; 5],

    solid_strings: Vec<(usize, usize)>,
    degenerate_letters: Vec<(usize, usize)>,
}

impl EDT {
    pub fn from_file(file_path: &str) -> Self {
        // TODO: read line by line?
        let data =
            fs::read(file_path).unwrap_or_else(|_| panic!("Unable to read file {}", file_path));
        Self::from_u8_slice(&data)
    }

    pub fn from_str(eds: &str) -> Self {
        Self::from_u8_slice(eds.as_bytes())
    }

    pub fn from_u8_slice(eds: &[u8]) -> Self {
        let col = Vec::<u8>::with_capacity(EXPECTED_ROWS);
        let mut data = vec![col; EXPECTED_COLS];

        let mut edges = Edges::new();

        // State variables
        let mut current_letter: Option<Letter> = None;
        let mut current_char: Option<Char> = None;

        let mut prev_letter: Option<Letter> = None;
        let mut prev_char: Option<Char> = None;

        let mut p = 0;
        let mut p_prime = 0;
        let mut max_p_prime = 0;
        let mut z = 0;
        let mut z_prime = 0;

        let mut length = 0;
        let mut size = 0;

        let mut solid_strings = Vec::<(usize, usize)>::with_capacity(EXPECTED_COLS);
        let mut degenerate_letters = Vec::<(usize, usize)>::with_capacity(EXPECTED_COLS);

        let mut solid_string_start: usize = 0 as usize;
        let mut solid_string_end: usize = 0;

        // A lookup for start indices of each nucleotide
        let mut start_indices = [
            HashSet::<usize>::new(), // A
            HashSet::<usize>::new(), // T
            HashSet::<usize>::new(), // C
            HashSet::<usize>::new(), // G
            HashSet::<usize>::new(), // epsilon/wildcard
        ];

        // Set the character and position of the nucleotide
        let mut set_lookup_index = |c: &u8, p: usize| {
            let mut lookup_index: usize = 0;
            match *c {
                b'A' => {}
                b'T' => lookup_index = 1,
                b'C' => lookup_index = 2,
                b'G' => lookup_index = 3,
                b'*' => lookup_index = 4,
                _ => panic!(),
            };

            start_indices[lookup_index].insert(p);
        };

        for c in eds {
            // automaton
            match c {
                b'A' | b'T' | b'C' | b'G' => {
                    prev_char = current_char;
                    prev_letter = current_letter;

                    current_char = Some(Char::Nucleotide);
                    if current_letter.is_none() {
                        // we started at a solid string
                        current_letter = Some(Letter::Solid);
                    }
                }

                b',' => {
                    // does not update the letter
                    prev_char = current_char;
                    current_char = Some(Char::Comma);
                }
                b'{' => {
                    prev_char = current_char;
                    prev_letter = current_letter;

                    current_char = Some(Char::Open);
                    current_letter = Some(Letter::Degenerate);
                }

                b'}' => {
                    prev_char = current_char;
                    prev_letter = current_letter;

                    current_char = Some(Char::Close);
                    current_letter = Some(Letter::Solid);
                }

                // handle whitespace
                b'\t' | b'\n' | b' ' => {
                    continue;
                }

                _ => panic!("Malformed EDS {}", *c as char),
            }

            // ---------------------------------
            // Helpers to query automaton state
            // ---------------------------------
            let is_solid_string_end = || -> bool {
                prev_char == Some(Char::Nucleotide) && current_char == Some(Char::Open)
            };

            // the start of a solid string
            let is_solid_string_start = || -> bool {
                // avoids the case of two degenerate letters
                (prev_char == Some(Char::Close) && current_char == Some(Char::Nucleotide)) ||
                // initial solid string
                    (data.len() == 0 &&
                     current_char == Some(Char::Nucleotide) &&
                     current_letter == Some(Letter::Solid))
            };

            // continuing a seed
            let in_degenerate_letter = || -> bool {
                current_letter == Some(Letter::Degenerate)
                    || current_char == Some(Char::Open)
                    || current_char == Some(Char::Close)
            };

            let in_solid_letter = || -> bool {
                current_letter == Some(Letter::Solid) && current_char == Some(Char::Nucleotide)
            };

            if is_solid_string_start() {
                solid_string_start = p;
            }

            if is_solid_string_end() {
                solid_string_end = p;
                solid_strings.push((solid_string_start, solid_string_end));
            }

            // ------------
            // Update State
            // ------------
            // solid_strings
            if in_solid_letter() {
                data[p] = vec![*c];
                set_lookup_index(c, p);
                p += 1;
                size += 1;
                length += 1;
            }

            // degenerate letters
            if in_degenerate_letter() {
                match (prev_char, current_char) {
                    // Case 1
                    // leading or ... empty degenerate letter
                    (Some(Char::Comma), Some(Char::Comma))
                    | (Some(Char::Open), Some(Char::Comma)) => {
                        // entire empty
                        p_prime = p;
                        z_prime += 1;
                    }

                    // Case 2
                    (Some(Char::Open), Some(Char::Nucleotide)) => {
                        p_prime = p;

                        match data.get_mut(p_prime) {
                            Some(col) => col.push(*c),
                            _ => data[p_prime] = vec![*c],
                        }

                        set_lookup_index(c, p_prime);

                        p_prime += 1;
                        size += 1;
                    }
                    // Case 3
                    (_, Some(Char::Comma)) | (_, Some(Char::Open)) => {
                        z_prime += 1;
                        p_prime = p;
                    }
                    // Case 4
                    (_, Some(Char::Close)) => {
                        // also handles trailing empty degenerate letter
                        // thanks to case 3

                        if PAD {
                            for col in p..p + max_p_prime {
                                for row in 0..z_prime {
                                    match data[col].get(row) {
                                        None => {
                                            data[col].push(WILDCARD);
                                            set_lookup_index(&b'*', col);
                                        }
                                        Some(_) => {}
                                    };
                                }
                            }
                        }

                        if z_prime > z {
                            z = z_prime;
                        }
                        degenerate_letters.push((p, p + max_p_prime));
                        z_prime = 0;
                        p = p + max_p_prime;
                        p_prime = 0;
                        max_p_prime = 0;
                        length += 1;

                        continue; // leave early
                    }
                    // Case 5
                    (_, Some(Char::Nucleotide)) => {
                        match data.get_mut(p_prime) {
                            Some(col) => col.push(*c),
                            _ => data[p_prime] = vec![*c],
                        }
                        set_lookup_index(c, p_prime);
                        p_prime += 1;
                        size += 1;
                    }
                    // case 6
                    _ => panic!("Unexpected case in degenerate letter"),
                }

                if p_prime - p > max_p_prime {
                    max_p_prime = p_prime - p
                }
            }
        }

        // add final offset info
        if current_char == Some(Char::Nucleotide) {
            solid_strings.push((solid_string_start, p));
        }

        data.truncate(p);
        degenerate_letters.shrink_to_fit();
        solid_strings.shrink_to_fit();
        data.shrink_to_fit();

        Self {
            data,
            edges,
            z,
            p,
            size,
            length,
            start_indices,
            degenerate_letters,
            solid_strings,
        }
    }

    pub fn print_stats(&self) -> String {
        format!(
            "z={} p={} length={} size={}\n\
             solid strings={:?}\n\
             degenerate letters={:?}\n\
             start indices={{\n\
             \tA={:?}\n\
             \tT={:?}\n\
             \tC={:?}\n\
             \tG={:?}\n\
             \t*={:?}\n\
             }}",
            self.z,
            self.p,
            self.length,
            self.size,
            self.solid_strings,
            self.degenerate_letters,
            self.start_indices[0],
            self.start_indices[1],
            self.start_indices[2],
            self.start_indices[3],
            self.start_indices[4],
        )
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn length(&self) -> usize {
        self.length
    }

    // should this be renamed diameter?
    pub fn p(&self) -> usize {
        self.p
    }

    // should this be renamed diameter?
    pub fn z(&self) -> usize {
        self.z
    }

    /// Positions containing a given char
    /// Allowed in lookup A, T, C, G, or * as u8
    pub fn get_start_indices(&self, c: u8) -> &HashSet<usize> {
        let mut lookup_index: usize = 0;
        match c {
            b'A' => {}
            b'T' => lookup_index = 1,
            b'C' => lookup_index = 2,
            b'G' => lookup_index = 3,
            _ => panic!("EDT::get_start_indices expects A, T, C or G. Got {c}"),
        };

        &self.start_indices[lookup_index]
    }

    pub fn get_solid_intervals(&self) -> &Vec<(usize, usize)> {
        &self.solid_strings
    }
    pub fn get_degenerate_letters(&self) -> &Vec<(usize, usize)> {
        &self.degenerate_letters
    }
}

impl Index<usize> for EDT {
    type Output = Vec<u8>;

    /// An index into the underlying vector
    ///
    /// ```text
    /// AT{TCC,C,}AA
    /// ```
    /// as
    ///
    /// ```text
    /// 0123456
    /// ATTCCAA
    ///   C**
    ///   ***
    /// ```

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx]
    }
}

impl Display for EDT {
    /// prints
    /// ```text
    /// A{CAT,GA,T}AT{TCC,C,}AA
    /// ```
    /// as
    /// ```text
    /// ACATATTCCAA
    ///  GA*  C**
    ///  T**  ***
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let rows = self.z;
        let cols = self.p;

        let mut s = String::new();
        for row in 0..rows {
            for col in 0..cols {
                match self.data[col].get(row) {
                    Some(r) => s.push(*r as char),
                    _ => s.push(' '),
                };
            }
            s.push('\n');
        }

        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    mod parse_str {
        use super::super::*;

        #[test]
        fn test_adjacent_degenerate_letters() {
            let ed_string = "A{T,G}{C,A}{T,A}TC";
            let edt = EDT::from_str(ed_string);
            assert_eq!(
                *edt.get_degenerate_letters(),
                Vec::from([(1, 2), (2, 3), (3, 4)])
            );

            let ed_string = "TGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GG";
            let edt = EDT::from_str(ed_string);
            assert_eq!(
                *edt.get_degenerate_letters(),
                Vec::from([(4, 5), (9, 10), (10, 11), (11, 12), (13, 14)])
            );
        }

        #[test]
        fn test_ignore_whitespace() {
            // spaces, tabs and newlines
            let ed_string = "ATCGATGGG{T,C}\n\
                             AACTT{TGTA,GGTA,C}AACTT{T, ,G}AG\t\
                             AG{G,T}CCGGTTTATATTGAT{T,C}CCTA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(
                *edt.get_degenerate_letters(),
                Vec::from([(9, 10), (15, 19), (24, 25), (29, 30), (45, 46)])
            );
        }

        #[test]
        fn test_degenerate_letters() {
            // -------
            // leading
            // -------
            let ed_string = "{TCGA,CTA,A}ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTT\
                             ATATTGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTG\
                             CTTGCTGTTG{A,G}CTC{T,G}TGAGTGAGCTTGCGAGATA";

            let edt = EDT::from_str(ed_string);
            assert_eq!(
                *edt.get_degenerate_letters(),
                Vec::from([
                    (0, 4),
                    (13, 14),
                    (19, 20),
                    (22, 23),
                    (38, 39),
                    (43, 44),
                    (44, 45),
                    (45, 46),
                    (47, 48),
                    (70, 71),
                    (74, 75)
                ])
            );

            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(0, 3), (5, 8)]);
            assert_eq!(*edt.get_degenerate_letters(), expected);

            // ------
            // middle
            // ------
            let ed_string = "AC{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(2, 5), (7, 10)]);
            assert_eq!(*edt.get_degenerate_letters(), expected);

            // --------
            // trailing
            // --------
            let ed_string = "AC{CAT,C,}AT{TCC,C,}AA{CGA,TC}";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(2, 5), (7, 10), (12, 15)]);
            assert_eq!(*edt.get_degenerate_letters(), expected);
        }

        #[test]
        fn test_empty_seed() {
            // trailing empty seed
            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(*edt.get_degenerate_letters(), Vec::from([(2, 5)]));

            // leading empty seed
            let ed_string = "AT{,TCC,C}AA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(*edt.get_degenerate_letters(), Vec::from([(2, 5)]));

            // TODO: what is this?
            // combined
            let ed_string = "AT{,,,}AA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(*edt.get_degenerate_letters(), Vec::from([(2, 2)]));
        }

        #[test]
        fn test_solid_strings() {
            // leading and trailing
            // multiple degenerate letters
            // starts with solid string
            let ed_string = "ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATATTGAT{T,C}CCTA\
                             {T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGTTG{A,G}CTC{T,G}\
                             TGAGTGAGCTTGCGAGATA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(
                *edt.get_solid_intervals(),
                Vec::from([
                    (0, 9),
                    (10, 15),
                    (16, 18),
                    (19, 34),
                    (35, 39),
                    (42, 43),
                    (44, 66),
                    (67, 70),
                    (71, 90)
                ])
            );

            // trailing
            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(3, 5), (8, 10)]);
            assert_eq!(*edt.get_solid_intervals(), expected);
        }
    }

    mod dag_construction {

        /*

        use super::super::*;

        #[ignore]
        #[test]
        fn test_edges() {
            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let out_edges = HashSet::from([6, 10, 9]);
            // dbg!(edt[5] as char, edt.edges);
            assert_eq!(*edt.from(5), out_edges);

            let in_edges = HashSet::from([5, 8, 9]);
            // assert_eq!(*edt.to(10), in_edges);

            let ed_string = "ATCGAAT{C,A}GAT{C,CATGC,,A}GA";
            let edt = EDT::from_str(ed_string);

            let out_edges = HashSet::from([7, 8]);
            assert_eq!(*edt.from(6), out_edges);

            let out_edges = HashSet::from([12, 13, 18, 19]);
            assert_eq!(*edt.from(11), out_edges);

            let in_edges = HashSet::from([7, 8]);
            assert_eq!(*edt.to(9), in_edges);

            let in_edges = HashSet::from([11, 12, 17, 18]);
            assert_eq!(*edt.to(19), in_edges);
        }
         */
    }

    // positions and offsets
    mod indexing {
        use super::super::*;

        // find char approximate char positions
        // echo "{CAT,C,}AT{TCC,C,}AA" | grep -aob '}'

        #[test]
        fn test_index_positions() {
            let ed_string = "A{CAT,GA,T}AT{TCC,C}AA";
            let edt = EDT::from_str(ed_string);

            let index = 2;
            let positions: Vec<char> = edt[index].iter().map(|e| *e as char).collect();
            assert_eq!(vec!['A', 'A', '*'], positions);

            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let index = 0;
            let positions: Vec<char> = edt[index].iter().map(|e| *e as char).collect();
            assert_eq!(positions, vec!['A']);

            let index = 4;
            let positions: Vec<char> = edt[index].iter().map(|e| *e as char).collect();
            assert_eq!(positions, vec!['C', '*', '*']);

            let index = edt.p() as usize - 1;
            let positions: Vec<char> = edt[index].iter().map(|e| *e as char).collect();
            assert_eq!(positions, vec!['A']);
        }

        #[test]
        fn test_find_start_indices() {
            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let start_a = HashSet::from([0, 6, 5]);
            assert_eq!(*edt.get_start_indices(b'A'), start_a);

            let start_g = HashSet::new();
            assert_eq!(*edt.get_start_indices(b'G'), start_g);
        }
    }

    mod matrix {
        use super::super::*;

        #[test]
        fn test_adjacent_degenerate_letters() {
            let ed_string = "A{T,G}{C,A}{T,A}TC";
            let edt = EDT::from_str(ed_string);
            eprintln!("{}", edt);
            eprintln!("{}\n", edt.print_stats());

            let ed_string = "TTATCTAC{CGCC,G}AC{CAAG,TC}G{C,G}AC{T,CGTA}{GTGT,G}\
                             TTCTC{GAAA,TG}ATATACC{AA,CG}GGTTCA{TAGC,GATT}CTTC\
                             {A,T}TAGGACTTCAGGGT{G,A}CAATT";
            let edt = EDT::from_str(ed_string);
            eprintln!("{}", edt);
            eprintln!("{}\n", edt.print_stats());

            // space
            let ed_string = "GA{T,C}ATT{T, ,G}ATC{,TGA,CA}GTA{GGCA,AGA,}CTT";
            let edt = EDT::from_str(ed_string);
            eprintln!("{}", edt);
            eprintln!("{}\n", edt.print_stats());

            let ed_string = "{TCGA,CTA,A}ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATAT\
                             TGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGT\
                             TG{A,G}CTC{T,G}TGAGTGAGCTTGCGAGATA";
            let edt = EDT::from_str(ed_string);
            eprintln!("{}", edt);
            eprintln!("{}\n", edt.print_stats());

            let ed_string = "{TCGA,CTA,A}ATCGATGGG{T,C}";
            let edt = EDT::from_str(ed_string);
            eprintln!("{}", edt);
            eprintln!("{}\n", edt.print_stats());
        }
    }
}

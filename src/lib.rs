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

 */

use std::collections::HashSet;
use std::fmt::{self, Display};
use std::fs;
use std::ops::Index;

// ---------
// Constants
// ----------

// Minimum number of chars it can hold without re-allocating
// expect an EDT to have a diameter of 50_000 bases
const EXPECTED_COLS: usize = 50_000;
const EXPECTED_ROWS: usize = 5; // expect a max of 5 variations
/// epsilon
pub const WILDCARD: u8 = b'*'; // epsilons
const A: u8 = b'A';
const T: u8 = b'T';
const C: u8 = b'C';
const G: u8 = b'G';
const PAD: bool = true; // pad the EDT with epsilons
const ENABLE_EDGES: bool = true;

// -------
// Aliases
// -------
/// A `[column, row]` index
pub type Coordinate = [usize; 2];

// ------------
// Helper types
// ------------
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

pub enum Edge {
    Incoming,
    Outgoing,
}

// ----------
// Traits
// ----------

pub trait Sequence {
    fn to_msa(&self) {}
}

// ----------
// Main types
// ----------

/// Contents of an [EDT](self::EDT) matrix cell
#[derive(Debug, Clone)]
pub struct Item {
    base: u8,
    incoming: HashSet<Coordinate>,
    outgoing: HashSet<Coordinate>,
}

impl Item {
    pub fn new(base: &u8) -> Self {
        Self {
            base: *base,
            outgoing: HashSet::<Coordinate>::new(),
            incoming: HashSet::<Coordinate>::new(),
        }
    }

    pub fn add_edge(&mut self, col: usize, row: usize, edge: Edge) -> Result<(), String> {
        match edge {
            Edge::Incoming => {
                self.incoming.insert([col, row]);
            }
            Edge::Outgoing => {
                self.outgoing.insert([col, row]);
            }
        }
        Ok(())
    }

    pub fn add_edges(&mut self, edges: HashSet<Coordinate>, edge: Edge) -> Result<(), String> {
        let successful = match edge {
            Edge::Incoming => edges
                .iter()
                .map(|e| self.incoming.insert(*e))
                .fold(true, |acc, v| acc && v),
            Edge::Outgoing => edges
                .iter()
                .map(|e| self.outgoing.insert(*e))
                .fold(true, |acc, v| acc && v),
        };

        if successful {
            Ok(())
        } else {
            Err(String::from("Fail"))
        }
    }

    pub fn base(&self) -> u8 {
        self.base
    }
}

/// Degenerate Text
pub struct DT {
    pub data: Vec<Vec<u8>>,
    pub z: usize,
}

impl DT {
    pub fn p(&self) -> usize {
        self.data.len()
    }

    pub fn z(&self) -> usize {
        self.z
    }
}

impl Index<usize> for DT {
    type Output = Vec<u8>;

    /// An index into the underlying vector
    ///
    /// ```text
    /// AT{T,C}AA
    /// ```
    /// as
    ///
    /// ```text
    /// 01234
    /// ATTAA
    ///   C
    /// ```

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx]
    }
}

impl Display for DT {
    /// prints
    /// ```text
    /// A{CAT,GAG}AT{TCC,CAA}AA
    /// ```
    /// as
    /// ```text
    /// ACATATTCCAA
    ///  GAG  CAA
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let cols = self.p();
        let rows = self.z();

        let mut s = String::new();
        for row in 0..rows {
            for col in 0..cols {
                match self.data[col].get(row) {
                    Some(i) => s.push(*i as char),
                    _ => s.push(' '),
                };
            }
            s.push('\n');
        }

        write!(f, "{}", s)
    }
}

impl Sequence for DT {
    fn to_msa(&self) {
        let rows = self.z;
        let cols = self.data.len();

        // let mut s = String::new();
        for (index, row) in (0..rows).enumerate() {
            println!(">{}", index);
            for col in 0..cols {
                match self.data[col].get(row) {
                    Some(i) => {
                        print!("{}", *i as char)
                    }
                    _ => {
                        print!("{}", *self.data[col].get(0).unwrap() as char)
                    }
                };
            }

            println!();
        }
    }
}

/// Elastic Degenerate Text
///
/// The underlying Elastic Degenerate Text
/// Can be thought of as a matrix
#[derive(Debug)]
pub struct EDT {
    /// the underlying matrix that holds the EDT
    data: Vec<Vec<Item>>,

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

    /// Intervals of solid strings
    solid_strings: Vec<(usize, usize)>,
    /// Intervals of degenerate letters
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
        let col = Vec::<Item>::with_capacity(EXPECTED_ROWS);
        let mut data = vec![col; EXPECTED_COLS];

        // State variables
        let mut current_letter: Option<Letter> = None;
        let mut current_char: Option<Char> = None;

        let mut prev_char: Option<Char>;

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
        let mut solid_string_end: usize;

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
                A => {}
                T => lookup_index = 1,
                C => lookup_index = 2,
                G => lookup_index = 3,
                WILDCARD => lookup_index = 4,
                _ => panic!(),
            };

            start_indices[lookup_index].insert(p);
        };

        for c in eds {
            // automaton
            match c {
                b'A' | b'T' | b'C' | b'G' => {
                    prev_char = current_char;

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

                    current_char = Some(Char::Open);
                    current_letter = Some(Letter::Degenerate);
                }

                b'}' => {
                    prev_char = current_char;

                    current_char = Some(Char::Close);
                    current_letter = Some(Letter::Solid);
                }

                // handle whitespace
                b'\t' | b'\n' | b' ' => {
                    continue;
                }

                _ => panic!("Malformed EDS {}", *c as char),
            }

            // -------------
            // Space related
            // -------------

            // reallocate the data vector
            if p >= data.len() - EXPECTED_COLS {
                let col = Vec::<Item>::with_capacity(EXPECTED_ROWS);
                let limit = p * 2;
                let data_slice = vec![col; limit];
                data.extend_from_slice(&data_slice);
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
            // Helpers
            // ------------
            let find_non_epsilon =
                |col: usize, row: usize, data: &Vec<Vec<Item>>| -> Option<Coordinate> {
                    // we get here if we have adj degenerate letters with epsilons
                    // go back and find first non wildcard
                    (0..col).rev().find_map(|col_index: usize| {
                        if data[col_index].len() == 1 {
                            // we get here if we have adj degenerate letters with completely empty variants
                            Some([col_index, 0])
                        } else if data[col_index][row].base() != WILDCARD {
                            Some([col_index, row])
                        } else {
                            None
                        }
                    })
                };

            let compute_incoming_edges = |p: usize, data: &Vec<Vec<Item>>| -> HashSet<Coordinate> {
                if p == 0 {
                    return HashSet::<Coordinate>::new();
                }

                let prev_col_index: usize = p - 1;
                let prev_col_len = data[prev_col_index].len();

                let incoming: HashSet<Coordinate> = (0..prev_col_len)
                    .filter_map(|row| {
                        if data[prev_col_index][row].base() != WILDCARD {
                            Some([prev_col_index, row])
                        } else {
                            find_non_epsilon(prev_col_index, row, &data)
                        }
                    })
                    .collect();

                incoming
            };

            // ------------
            // Update State
            // ------------
            // solid_strings
            if in_solid_letter() {
                let mut item = Item::new(c);

                if ENABLE_EDGES {
                    let incoming: HashSet<Coordinate> = compute_incoming_edges(p, &data);
                    item.add_edges(incoming.clone(), Edge::Incoming)
                        .expect(&format!(
                            "Failed in solid_string to add incoming edges p={} s={}",
                            p, size
                        ));
                    for [col, row] in incoming {
                        // No need to check for wildcard as it was checked by incoming
                        data[col][row]
                            .add_edge(p, 0, Edge::Outgoing)
                            .expect(&format!(
                                "Failed in solid_string to add outgoing edges p={} s={}",
                                p, size
                            ));
                    }
                }
                data[p] = vec![item];
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
                        let mut item = Item::new(c);

                        if ENABLE_EDGES {
                            let incoming: HashSet<Coordinate> =
                                compute_incoming_edges(p_prime, &data);
                            item.add_edges(incoming.clone(), Edge::Incoming)
                                .expect(&format!(
                                    "Failed case 2 to add incoming edges p={} s={}",
                                    p_prime, size
                                ));
                            for [col, row] in incoming {
                                // No need to check for wildcard as it was checked by incoming
                                data[col][row]
                                    .add_edge(p_prime, z_prime, Edge::Outgoing)
                                    .expect(&format!(
                                        "Failed case 2 to add outgoing edges p={} s={}",
                                        p_prime, size
                                    ));
                            }
                        }

                        match data.get_mut(p_prime) {
                            Some(col) => col.push(item),
                            _ => data[p_prime] = vec![item],
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
                                            data[col].push(Item::new(&WILDCARD));
                                            set_lookup_index(&WILDCARD, col);
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
                        let mut item = Item::new(c);

                        if ENABLE_EDGES {
                            let incoming: HashSet<Coordinate> =
                                compute_incoming_edges(p_prime, &data);
                            item.add_edges(incoming.clone(), Edge::Incoming)
                                .expect(&format!(
                                    "Failed case 5 to add incoming edges p={} s={}",
                                    p_prime, size
                                ));

                            for [col, row] in incoming {
                                // No need to check for wildcard as it was checked by incoming
                                data[col][row]
                                    .add_edge(p_prime, z_prime, Edge::Outgoing)
                                    .expect(&format!(
                                        "Failed case 5 to add outgoing edges p={} s={}",
                                        p_prime, size
                                    ));
                            }
                        }

                        match data.get_mut(p_prime) {
                            Some(col) => col.push(item),
                            _ => data[p_prime] = vec![item],
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

    pub fn base_at(&self, index: Coordinate) -> u8 {
        let [col, row] = index;
        self.data[col][row].base()
    }

    pub fn outgoing(&self, index: Coordinate) -> &HashSet<Coordinate> {
        let [col, row] = index;
        &self.data[col][row].outgoing
    }

    pub fn incoming(&self, index: Coordinate) -> &HashSet<Coordinate> {
        let [col, row] = index;
        &self.data[col][row].incoming
    }

    /// Positions containing a given char
    /// Allowed in lookup A, T, C, G, or * as u8
    pub fn get_start_indices(&self, c: u8) -> &HashSet<usize> {
        let mut lookup_index: usize = 0;
        match c {
            A => {}
            T => lookup_index = 1,
            C => lookup_index = 2,
            G => lookup_index = 3,
            WILDCARD => lookup_index = 3,
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

    pub fn extract_inelastic(&self) -> DT {
        let diagonal = self.p();
        let mut z = 0;
        let mut degenerate_matrix: Vec<Vec<u8>> = Vec::with_capacity(diagonal);

        for index in 0..diagonal {
            let col: &Vec<Item> = &self[index];
            let bases_in_col: Vec<u8> = col.iter().map(|item: &Item| item.base()).collect();
            if !bases_in_col.iter().copied().any(|c: u8| c == WILDCARD) {
                let z_prime = bases_in_col.len();
                if z_prime > z {
                    z = z_prime;
                }
                degenerate_matrix.push(bases_in_col);
            }
        }

        DT {
            data: degenerate_matrix,
            z,
        }
    }
}

impl Sequence for EDT {
    fn to_msa(&self) {
        let rows = self.z;
        let cols = self.p;

        // let mut s = String::new();
        for (index, row) in (0..rows).enumerate() {
            println!(">{}", index);
            for col in 0..cols {
                match self.data[col].get(row) {
                    Some(i) => {
                        let c = if WILDCARD == i.base() {
                            '-'
                        } else {
                            i.base() as char
                        };
                        print!("{}", c)
                    }
                    _ => {
                        print!("{}", self.data[col].get(0).unwrap().base() as char)
                    }
                };
            }

            println!();
        }

        // write!(f, "{}", s)
    }
}

impl Index<usize> for EDT {
    type Output = Vec<Item>;

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
                    Some(i) => s.push(i.base() as char),
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
            let positions: Vec<char> = edt[index].iter().map(|e| e.base() as char).collect();
            assert_eq!(vec!['A', 'A', '*'], positions);

            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let index = 0;
            let positions: Vec<char> = edt[index].iter().map(|e| e.base() as char).collect();
            assert_eq!(positions, vec!['A']);

            let index = 4;
            let positions: Vec<char> = edt[index].iter().map(|e| e.base() as char).collect();
            assert_eq!(positions, vec!['C', '*', '*']);

            let index = edt.p() as usize - 1;
            let positions: Vec<char> = edt[index].iter().map(|e| e.base() as char).collect();
            assert_eq!(positions, vec!['A']);
        }

        #[test]
        fn test_find_start_indices() {
            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let start_a = HashSet::from([0, 6, 5]);
            assert_eq!(*edt.get_start_indices(A), start_a);

            let start_g = HashSet::new();
            assert_eq!(*edt.get_start_indices(G), start_g);
        }
    }

    mod dag {
        use super::super::*;

        #[test]
        fn test_incoming_edges() {
            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            if ENABLE_EDGES {
                // solid
                let in_edges = HashSet::from([[8, 0]]);
                assert_eq!(*edt.incoming([9, 0]), in_edges);

                let in_edges = HashSet::from([[7, 0], [5, 1], [4, 0]]);
                assert_eq!(*edt.incoming([8, 0]), in_edges);

                let in_edges = HashSet::from([[4, 0]]);
                assert_eq!(*edt.incoming([5, 0]), in_edges);

                // degenerate
                let in_edges = HashSet::from([[4, 0]]);
                assert_eq!(*edt.incoming([5, 1]), in_edges);

                let in_edges = HashSet::from([]);
                assert_eq!(*edt.incoming([5, 2]), in_edges);

                let in_edges = HashSet::from([[6, 0]]);
                assert_eq!(*edt.incoming([7, 0]), in_edges);
            }
        }

        #[test]
        fn test_outgoing_edges() {
            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            if ENABLE_EDGES {
                // solid
                let out_edges = HashSet::<Coordinate>::new();
                assert_eq!(*edt.outgoing([9, 0]), out_edges);

                let out_edges = HashSet::from([[9, 0]]);
                assert_eq!(*edt.outgoing([8, 0]), out_edges);

                // degenerate
                let out_edges = HashSet::from([[3, 0]]);
                assert_eq!(*edt.outgoing([0, 1]), out_edges);

                let out_edges = HashSet::from([[8, 0]]);
                assert_eq!(*edt.outgoing([5, 1]), out_edges);
            }
        }
    }

    mod dt {
        use super::super::*;
        #[test]
        fn test_foo() {
            let ed_string = "A{T,G}{C,A}{T,A}TC";
            let edt = EDT::from_str(ed_string);

            let dt = edt.extract_inelastic();

            dt.to_msa();
        }
    }
}

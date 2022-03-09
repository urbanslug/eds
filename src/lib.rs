/*!
EDS parser

Parse an EDS into a DAG of its characters.

Degenerate letters are between brackets.
Individual seeds are separated by commas.
Empty seeds are allowed
`AT{TCC,AG,C,}AA`

Generate EDS from fasta file and VCF
[https://github.com/webmasterar/edso]()

Generate random EDS
[https://github.com/webmasterar/EDSRand]()

```
use eds::EDT;
use std::collections::HashSet;

let ed_string = "AT{TCC,C,}AA";
let edt = EDT::from_str(ed_string);

// size and length
assert_eq!(edt.size() as usize, 8);
assert_eq!(edt.length() as usize, 3);

// edges
let out_edges = HashSet::from([2, 5, 6]);
assert_eq!(*edt.from(1), out_edges);

// indexing
assert_eq!(edt[0], b'A');
assert_eq!(edt[5], b'C');
assert_eq!(edt[( edt.size() as usize - 1) ], b'A');
```
 */

use std::ops::Index;
use std::collections::HashSet;
use std::iter::IntoIterator;

// Minimum number of chars it can hold without re-allocating
// assume small viral pan-genome (50k)
const EXPECTED_MIN_SIZE: usize = 50_000;

/// 0 => from (incoming nodes)
/// 1 => to (outgoing nodes)
#[derive(Debug)]
struct Edges(HashSet<u32>, HashSet<u32>);

impl Edges {
    fn new() -> Self {
        Edges(HashSet::<u32>::new(), HashSet::<u32>::new())
    }

    // outgoing edges
    fn from(&self) -> &HashSet<u32> {
        &self.1
    }

    // incoming edges
    fn to(&self) -> &HashSet<u32> {
        &self.0
    }
}

impl<'a> Iterator for &'a Edges {
    type Item = &'a HashSet<u32>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(&self.1)
    }
}
/// A letter in the EDT. It can either be degenerate or solid
#[derive(Debug, PartialEq, Copy,  Clone)]
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
#[derive(Debug, PartialEq, Copy,  Clone)]
enum Char {
    Nucleotide,
    Comma,
    Open,
    Close
}

#[derive(Debug)]
pub struct EDT {
    /// The chars of the EDT
    data: Vec<u8>,

    /// Connections between chars
    // TODO: use an adj matrix
    // In edges ... out edges
    edges: Vec<Edges>,

    /// size (N)
    length: u32,

    /// length (n)
    size: u32,

    /// An array of sets of start indices
    /// Start indices for each nucleotide
    /// Index of where each char starts
    /// A = 0 T = 1 C = 2 G = 3
    start_indices: [HashSet<u32>; 4],

    solid_strings: Vec<(usize, usize)>,
    degenerate_letters: Vec<(usize, usize)>,
}

// ----------
// Helpers
// ----------
fn is_valid_index<T>(idx: usize, vec: &Vec<T>) -> bool {
    idx < vec.len()
}

impl EDT {
    // ----------
    // Methods
    // ----------
    #[allow(unused_assignments)]
    pub fn from_str(eds: &str) -> Self {
        let mut data = Vec::<u8>::with_capacity(EXPECTED_MIN_SIZE);
        let mut edges = Vec::<Edges>::with_capacity(EXPECTED_MIN_SIZE);
        let mut solid_strings = Vec::<(usize, usize)>::with_capacity(EXPECTED_MIN_SIZE);
        let mut degenerate_letters = Vec::<(usize, usize)>::with_capacity(EXPECTED_MIN_SIZE);
        // Lookup for start indices of each nucleotide
        let mut start_indices = [
            HashSet::<u32>::new(), // A
            HashSet::<u32>::new(), // T
            HashSet::<u32>::new(), // C
            HashSet::<u32>::new(), // G
        ];

        // set the character and position of the nucleotide
        let mut set_lookup_index = |c: &u8, size: u32| {
            let mut lookup_index: usize = 0;
            match *c {
                b'A' => {},
                b'T' => lookup_index = 1,
                b'C' => lookup_index = 2,
                b'G' => lookup_index = 3,
                _ => panic!()
            };

            start_indices[lookup_index].insert(size);
        };

        // the index of the current character
        let mut size: u32 = 0;

        // number of letters
        let mut length: u32 = 1;

        // State variables
        let mut current_letter: Option<Letter> = None;
        let mut current_char: Option<Char> = None;

        let mut prev_letter: Option<Letter> = None;
        let mut prev_char: Option<Char> = None;

        let mut seed_starts: Option<HashSet<u32>> = None;
        let mut seed_stops: Option<HashSet<u32>> = None;
        let mut solid_end: Option<u32> = None;
        let mut solid_start: Option<u32> = None;

        let mut solid_string_start: Option<u32> = None;
        let mut solid_string_end: Option<u32> = None;
        let mut degenerate_letter_start: Option<u32> = None;
        let mut degenerate_letter_end: Option<u32> = None;

        let mut solid_string_count: u32 = 0;
        let mut degenerate_letter_count: u32 = 0;
        let mut contains_empty_seed = false;

        for c in eds.as_bytes() {
            // automaton
            match c {
                b'A' | b'T'|  b'C' | b'G' => {
                    prev_char = current_char;
                    prev_letter = current_letter;

                    current_char = Some(Char::Nucleotide);
                    if current_letter.is_none() {
                        // we started at a solid string
                        current_letter = Some(Letter::Solid);
                    }
                },

                b',' => {
                    // does not update the letter
                    prev_char = current_char;
                    current_char = Some(Char::Comma);
                },
                b'{' => {
                    prev_char = current_char;
                    prev_letter = current_letter;

                    current_char = Some(Char::Open);
                    current_letter = Some(Letter::Degenerate);
                },

                b'}' => {
                    prev_char = current_char;
                    prev_letter = current_letter;

                    current_char = Some(Char::Close);
                    current_letter = Some(Letter::Solid);
                }

                _ => panic!("Malformed EDS {}", *c as char)
            }

            // ---
            // Inside a degenerate letter
            // ---
            let is_seed_start = || -> bool {
                prev_char == Some(Char::Open) ||
                    prev_char == Some(Char::Comma)
            };

            let is_seed_stop = || -> bool {
                current_char == Some(Char::Close) ||
                    prev_char == Some(Char::Comma)
            };

            let found_empty_seed = || -> bool {
                ( prev_char == Some(Char::Open) &&
                  current_char == Some(Char::Comma)
                ) || ( prev_char == Some(Char::Comma) &&
                        current_char == Some(Char::Comma) ) ||
                    ( prev_char == Some(Char::Comma) &&
                      current_char == Some(Char::Close)
                    )
            };

            // continuing a seed
            let cont_seed = || -> bool {
                current_letter == Some(Letter::Degenerate) &&
                    prev_letter == Some(Letter::Degenerate) &&
                    prev_char == Some(Char::Nucleotide) &&
                    current_char == Some(Char::Nucleotide)
            };

            let is_degenerate_letter_start = || -> bool {
                current_char == Some(Char::Open)
            };

            let is_degenerate_letter_end = || -> bool {
                current_char == Some(Char::Close)
            };

            // ---
            // Solid string
            // ---

            // from the perspective of the solid string
            // the end of a solid string
            // that is *not* at the very end of an EDS
            let is_solid_string_end = || -> bool {
                prev_char == Some(Char::Nucleotide) &&
                    current_char == Some(Char::Open)
            };

            // the start of a solid string
            let is_solid_string_start = || -> bool {
                current_char == Some(Char::Close) ||
                    (data.len() == 0 &&
                     current_char == Some(Char::Nucleotide) &&
                     current_letter == Some(Letter::Solid))
            };

            let is_first_char_in_leading_solid_string = || -> bool {
                data.len() == 0 &&
                 current_char == Some(Char::Nucleotide) &&
                 current_letter == Some(Letter::Solid)
            };

            // continuing a solid string
            let cont_solid = || -> bool {
                current_letter == Some(Letter::Solid) &&
                    prev_letter == Some(Letter::Solid) &&
                    prev_char == Some(Char::Nucleotide)
            };


            // ----
            // Update DAG state
            // ----

            if is_seed_start() {
                match seed_starts.as_mut() {
                    Some(s) => { s.insert(size); },
                    _ => {
                        let mut s = HashSet::<u32>::new();
                        s.insert(size);
                        seed_starts = Some(s);
                    }
                }
            }


            if is_seed_stop() {
                let seed_stop_idx = size - 1;
                match seed_stops.as_mut() {
                    Some(s) => { s.insert(seed_stop_idx); },
                    _ => {
                        let mut s = HashSet::<u32>::new();
                        s.insert(seed_stop_idx);
                        seed_stops = Some(s);
                    }
                }
            }

            if is_degenerate_letter_start() {
                degenerate_letter_start = Some(size);
                degenerate_letter_count += 1;
            }

            if is_degenerate_letter_end() {
                degenerate_letter_end = Some(size);
            }

            /*
            if is_solid_end() {
                match solid_end.as_mut() {
                    None => { solid_end = Some(size) },
                    _ => {},
                };
            }
            */

            if is_solid_string_end() {
                solid_string_end = Some(size-1);
                solid_end = Some(size-1);
            }

            if is_solid_string_start() {
                solid_string_start = Some(size);
                solid_string_count += 1;
            }

            if found_empty_seed() {
                contains_empty_seed = true;
            }

            // TODO: improve this
            let mut update_offsets = || {
                if is_degenerate_letter_end() {
                    degenerate_letters.push(
                        (degenerate_letter_start.unwrap() as usize, degenerate_letter_end.unwrap() as usize)
                    );

                    degenerate_letter_start = None;
                    degenerate_letter_end = None;
                }

                if is_solid_string_end() {
                    if solid_string_start.is_none() {
                        solid_strings.push((0, solid_string_end.unwrap() as usize + 1) )
                    } else {
                        solid_strings.push(
                            (solid_string_start.unwrap() as usize, solid_string_end.unwrap() as usize + 1)
                        )
                    }

                    solid_string_start = None;
                    solid_string_end = None;
                }
            };

            // update edges
            let mut update_edges = || {
                // no need to check that the current char is a nucleotide
                // we may have just enetered a solid string from a Char::Closing
                if cont_solid() || cont_seed()
                {
                    let curr: usize = size as usize;
                    let prev: usize = curr - 1;

                    // if it's not large enough just pad it
                    while !is_valid_index(curr, &edges) {
                        edges.push(Edges::new());
                    }

                    // set prev outgoing edges
                    edges[prev].1.insert(size);

                    // set the current incoming edges
                    edges[curr].0.insert(size-1);
                }

                if is_solid_string_start() {

                    // if it's not large enough just pad it
                    while !is_valid_index(size as usize, &edges) {
                        edges.push(Edges::new());
                    }

                    // add an edge between the two adj solid strings
                    if contains_empty_seed && solid_string_count > 1 && is_solid_string_start()  {
                        let from = solid_end.expect(&format!("{size}"));
                        let to_idx = size as usize;
                        let from_idx = from as usize;
                        edges[from_idx].1.insert(size);
                        edges[to_idx].0.insert(from);
                    }

                    if solid_string_count > 1 && solid_end.is_some() {
                        // add edge btwn seeds starts and previous solid strings
                        match seed_starts.as_ref() {
                            Some(s) => {
                                for i in s {
                                    edges[*i as usize].0.insert(solid_end.expect(&format!("{size}")));
                                    edges[solid_end.unwrap() as usize].1.insert(*i);
                                }
                            },
                            _ => {}
                        }
                    }

                    // add edge btwn seed ends and next solid strings
                    match seed_stops.as_ref() {
                        Some(s) => {
                            for i in s {
                                edges[*i as usize].1.insert(size);
                                edges[size as usize].0.insert(*i);
                            }
                        },
                        _ => {}
                    }


                    seed_starts = None;
                    seed_stops = None;
                    solid_end = None;
                    solid_start = None;
                }
            };

            // update data
            update_edges();
            update_offsets();

            if current_char == Some(Char::Nucleotide) {
                set_lookup_index(c, size);
                // save the current letter
                data.push(*c);
                // inc size
                size += 1;
            }

            // update letter
            // If we have changed the letter
            // and we were in a letter
            if prev_letter.is_some() && prev_letter != current_letter {
                length += 1;
            }

        }

        // add final offset info
        if current_char == Some(Char::Nucleotide) {
            solid_strings.push (
                (solid_string_start.unwrap() as usize, size as usize))
        }

        degenerate_letters.shrink_to_fit();
        solid_strings.shrink_to_fit();
        data.shrink_to_fit();
        edges.shrink_to_fit();

        EDT {
            data,
            edges,
            size,
            length,
            start_indices,
            degenerate_letters,
            solid_strings,
        }
    }

    /// Get the char at position
    pub fn get(&self, idx: usize) -> Option<u8> {
        self.data.get(idx).cloned()
    }

    pub fn get_data(&self) -> &Vec<u8> {
        &self.data
    }

    pub fn get_solid_string_offsets(&self) -> &Vec<(usize, usize)> {
        &self.solid_strings
    }

    pub fn get_degenerate_letter_offsets(&self) -> &Vec<(usize, usize)> {
        &self.degenerate_letters
    }

    /// Positions containing a given char
    /// Allowed in lookup A, T, C or G as u8
    pub fn get_start_indices(&self, c: u8) -> &HashSet<u32> {
        let mut lookup_index: usize = 0;
        match c {
            b'A' => {},
            b'T' => lookup_index = 1,
            b'C' => lookup_index = 2,
            b'G' => lookup_index = 3,
            _ => panic!("EDT::get_start_indices expects A, T, C or G. Got {c}")
        };

        &self.start_indices[lookup_index]
    }

    /// The outgoing edges from nucleotide at index idx
    pub fn from(&self, idx: usize) -> &HashSet<u32> {
        &self.edges[idx].from()
    }

    /// indices of the nodes with edges coming into nucleotide at index idx
    pub fn to(&self, idx: usize) -> &HashSet<u32> {
        &self.edges[idx].to()
    }

    pub fn size(&self) -> u32 {
        self.size
    }

    pub fn length(&self) -> u32 {
        self.length
    }
}

impl Index<usize> for EDT {
    type Output = u8;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx]
    }
}

// ---------
// Iterator
// ---------
impl IntoIterator for EDT {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

#[cfg(test)]
mod tests {
    mod parse_str {
        use super::super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn test_adjacent_degenerate_letters() {

            let ed_string = "CCTA{T,G}{T,A}CAG";
            let edt = EDT::from_str(ed_string);
            assert!(false);
            // assert_eq!(edt.size as usize, edt.data.len());

            let ed_string = "TGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GG";

            let edt = EDT::from_str(ed_string);
            assert!(false);
            // assert_eq!(edt.size as usize, edt.data.len());
        }

        #[test]
        fn test_start_degenerate() {
            let ed_string = "{TCGA,CTA,A}ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATAT\
                             TGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGT\
                             TG{A,G}CTC{T,G}TGAGTGAGCTTGCGAGATA";

            let edt = EDT::from_str(ed_string);
            assert_eq!(edt.size as usize, edt.data.len());
        }

        #[test]
        fn test_empty_seed() {
            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(edt.size as usize, edt.data.len());

            let ed_string = "AT{,TCC,C}AA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(edt.size as usize, edt.data.len());

            let ed_string = "AT{,,,}AA";
            let edt = EDT::from_str(ed_string);
            assert_eq!(edt.size as usize, edt.data.len());
        }

        #[test]
        fn test_start_solid() {

            // starts with solid string
            let ed_string = "AT{TCC,AG,C}AA";
            let edt = EDT::from_str(ed_string);

            assert_eq!(edt.size as usize, edt.data.len());

            // multiple degenerate letters
            // starts with solid string
            let ed_string = "ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATATTGAT{T,C}CCTA\
                             {T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGTTG{A,G}CTC{T,G}\
                             TGAGTGAGCTTGCGAGATA";
            let edt = EDT::from_str(ed_string);

            assert_eq!(edt.size as usize, edt.data.len());
        }

        #[test]
        fn test_multiple_degenerate_letters() {
            let ed_string = "ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATATTGAT{T,C}CCTA\
                             {T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGTTG{A,G}CTC{T,G}\
                             TGAGTGAGCTTGCGAGATA";
            let edt = EDT::from_str(ed_string);

            assert_eq!(edt.size as usize, edt.data.len());
        }
    }
    mod dag_construction {
        use super::super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn test_edges() {
            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let out_edges = HashSet::from([6, 10, 9]);
            // dbg!(edt[5] as char, edt.edges);
            assert_eq!(*edt.from(5), out_edges);

            let in_edges = HashSet::from([5,8,9]);
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

        #[test]
        fn test_indexing() {
            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            assert_eq!(edt[0], b'A');
            assert_eq!(edt[5], b'C');
            assert_eq!(edt[( edt.size() as usize - 1) ], b'A');
        }

        #[test]
        fn test_find_start_indices() {
            let ed_string = "AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);


            let start_a = HashSet::from([0, 6, 7]);
            assert_eq!(*edt.get_start_indices(b'A'), start_a);

            let start_g = HashSet::new();
            assert_eq!(*edt.get_start_indices(b'G'), start_g);
        }
        }
    mod offsets {
        use super::super::*;
        use pretty_assertions::assert_eq;

        // find char approximate char positions
        // echo "{CAT,C,}AT{TCC,C,}AA" | grep -aob '}'

        #[test]
        fn test_degenerate_letters() {
            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(0, 4), (6, 10)]);
            assert_eq!(edt.degenerate_letters, expected);

            let ed_string = "AC{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(2, 6), (8, 12)]);
            assert_eq!(edt.degenerate_letters, expected);

            let ed_string = "AC{CAT,C,}AT{TCC,C,}AA{CGA,TC}";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(2, 6), (8, 12), (14, 19)]);
            assert_eq!(edt.degenerate_letters, expected);
        }

        #[test]
        fn test_solid_strings() {


            let ed_string = "{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(4, 6), (10,12)]);
            assert_eq!(edt.solid_strings, expected);

            let ed_string = "AC{CAT,C,}AT{TCC,C,}AA";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(0,2), (6,8), (12, 14)]);
            assert_eq!(edt.solid_strings, expected);


            let ed_string = "AC{CAT,C,}AT{TCC,C,}AA{CGA,TC}";
            let edt = EDT::from_str(ed_string);

            let expected = Vec::from([(0,2), (6, 8), (12, 14)]);
            assert_eq!(edt.solid_strings, expected);


        }

    }
    mod iterators {
        use super::super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn test_adjacent_degenerate_letters() {
            let ed_string = "TGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GG";
            let _edt = EDT::from_str(ed_string);
        }
    }
}

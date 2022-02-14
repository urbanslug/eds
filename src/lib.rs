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
 */

use std::collections::HashSet;

/// 0 => from (incoming nodes)
/// 1 => to (outgoing nodes)
#[derive(Debug)]
struct Edges(HashSet<u32>, HashSet<u32>);

impl Edges {
    fn new() -> Self {
        Edges(HashSet::<u32>::new(), HashSet::<u32>::new())
    }

    fn from(&self) -> &HashSet<u32> {
        &self.0
    }

    fn to(&self) -> &HashSet<u32> {
        &self.0
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
}


impl EDT {
    // ----------
    // Helpers
    // ----------
    fn is_valid_index<T>(idx: usize, vec: &Vec<T>) -> bool {
        idx < vec.len()
    }

    // ----------
    // Methods
    // ----------
    #[allow(unused_assignments)]
    pub fn from_str(eds: &str) -> Self {

        let mut data = Vec::<u8>::new();
        // let empty_set = HashSet::<u32>::new();
        let mut edges: Vec<Edges> = Vec::new();

        // Lookup for start indices of each nucleotide
        let mut start_indices = [
            HashSet::<u32>::new(), // A
            HashSet::<u32>::new(), // T
            HashSet::<u32>::new(), // C
            HashSet::<u32>::new(), // G
        ];
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

        // Assume that we start in a solid string
        let mut length: u32 = 1;

        let mut current_letter: Option<Letter> = None;
        let mut current_char: Option<Char> = None;

        let mut prev_letter: Option<Letter> = None;
        let mut prev_char: Option<Char> = None;

        let mut seed_starts: Option<HashSet<u32>> = None;
        let mut seed_stops: Option<HashSet<u32>> = None;
        let mut solid_end: Option<u32> = None;
        let mut solid_start: Option<u32> = None;

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


            let is_seed_start = || -> bool {
                prev_char == Some(Char::Open) ||
                    prev_char == Some(Char::Comma)
            };

            let is_seed_stop = || -> bool {
                current_char == Some(Char::Close) ||
                    prev_char == Some(Char::Comma)
            };

            // the end of a solid string
            // that is *not* at the very end of an EDS
            let is_solid_end = || -> bool {
                current_char == Some(Char::Open) ||
                    prev_char == Some(Char::Nucleotide)
            };

            let has_empty_seed = || -> bool {
                (prev_char == Some(Char::Comma) &&
                 current_char == Some(Char::Comma) ) ||
                    ( prev_char == Some(Char::Comma) &&
                      current_char == Some(Char::Close)
                    ) ||
                    ( prev_char == Some(Char::Open) &&
                      current_char == Some(Char::Comma)
                    )
            };

            // the end of a solid string
            // that is *not* at the start of an EDS
            let is_solid_start = || -> bool {
                current_char == Some(Char::Close)
            };

            // continuing a seed
            let cont_seed = || -> bool {
                current_letter == Some(Letter::Degenerate) &&
                    prev_letter == Some(Letter::Degenerate) &&
                    prev_char == Some(Char::Nucleotide) &&
                    current_char == Some(Char::Nucleotide)
            };

            // continuing a solid string
            let cont_solid = || -> bool {
                current_letter == Some(Letter::Solid) &&
                    prev_letter == Some(Letter::Solid) &&
                    prev_char == Some(Char::Nucleotide)
            };


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

            if is_solid_end() {
                match solid_end.as_mut() {
                    None => { solid_end = Some(size) },
                    _ => {},
                }
            }


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

                if is_solid_start() {

                    // if it's not large enough just pad it
                    while !is_valid_index(size as usize, &edges) {
                        edges.push(Edges::new());
                    }

                    if has_empty_seed() {
                        let from = solid_end.unwrap();
                        let to_idx = size as usize;
                        let from_idx = from as usize;
                        edges[from_idx].1.insert(size);
                        edges[to_idx].0.insert(from);
                    }

                    match seed_starts.as_ref() {
                        Some(s) => {
                            for i in s {
                                edges[*i as usize].0.insert(solid_end.unwrap());
                                edges[solid_end.unwrap() as usize].1.insert(*i);
                            }
                        },
                        _ => {}
                    }

                    match seed_stops.as_ref() {
                        Some(s) => {
                            for i in s {
                                edges[*i as usize].1.insert(size);
                                edges[size as usize].0.insert(*i);
                            }
                        },
                        _ => {}
                    }


                    // reset
                    seed_starts = None;
                    seed_stops = None;
                    solid_end = None;
                    solid_start = None;
                }
            };

            update_edges();

            // update data
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

        EDT {
            data,
            edges,
            size,
            length,
            start_indices
        }
    }

    /// Get the char at position
    pub fn at(&self, idx: usize) -> Option<u8> {
        self.data.get(idx).cloned()
    }

    /// Positions containing a given char
    /// Allowed in lookup A, T, C or G as u8
    pub fn get_start_indices(&self, c: u8) -> Option<&HashSet<u32>> {
        let mut lookup_index: usize = 0;
        match c {
            b'A' => {},
            b'T' => lookup_index = 1,
            b'C' => lookup_index = 2,
            b'G' => lookup_index = 3,
            _ => { return None }
        };
        self.start_indices.get(lookup_index)
    }

    fn get_edges(&self, idx: usize) -> Option<&Edges> {
        self.edges.get(idx)
    }

    /// get the incoming edges
    pub fn from(&self, idx: usize) -> Option<&HashSet<u32>> {
        self.get_edges(idx).map(|x: &Edges| x.from())
    }

    /// get the outgoing edges
    pub fn to(&self, idx: usize) -> Option<&HashSet<u32>> {
        self.get_edges(idx).map(|x: &Edges| x.to())
    }

    pub fn size(&self) -> u32 {
        self.size
    }

    pub fn length(&self) -> u32 {
        self.length
    }
}

#[cfg(test)]
mod tests {
    mod parse_str {
        use super::super::*;

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
            // let ed_string = "AT{TCC,C}AA";
            let ed_string = "AT{TCC,AG,C}AA";
            // let ed_string = "AT{TCC,C,}AA";
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
}

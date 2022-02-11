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

#[derive(Debug ,PartialEq, Copy,  Clone)]
enum LetterType {
    SolidString,
    DegenerateLetter,
    SeedBoundary,
    StartElastic,
    EndElastic,
}

use LetterType::{*};

const EMPTY_EDGE_LIST: Vec<u32> = Vec::new();
const EMPTY_EDGE_PIAR: (Vec<u32>, Vec<u32>) = (EMPTY_EDGE_LIST, EMPTY_EDGE_LIST);

/// 0 vec of in edges
/// 1 vec of out edges
type EdgesPair = (Vec<u32>, Vec<u32>);

#[derive(Debug)]
pub struct EDT {
    /// The chars of the EDT
    data: Vec<u8>,

    /// Connections between chars
    // TODO: use an adj matrix
    // In edges ... out edges
    edges: Vec<EdgesPair>,

    /// size (N)
    length: u32,

    /// length (n)
    size: u32
}

impl EDT {
    pub fn from_str(eds: &str) -> Self {
        let mut data = Vec::<u8>::new();
        let mut edges: Vec<(Vec<u32>, Vec<u32>)> = Vec::new();

        // the index of the current character
        let mut counter: i32 = 0;

        // Assume that we start in a solid string
        let mut length: u32 = 1;
        let mut current_index_letter_type = SolidString;
        let mut prev_index_letter_type = SolidString;

        let mut last_solid_string_index: Option<u32> = None;
        let mut next_solid_string_index: Option<u32> = None;
        let mut seed_indices_start: Option<Vec<u32>> = None;
        let mut seed_indices_stop: Option<Vec<u32>> = None;

        for c in eds.as_bytes() {

            eprintln!(" {}\t{}\t{:?}\t\t{:?}", *c as char, counter, prev_index_letter_type, current_index_letter_type);

            let prev: i32 = counter - 1;

            match c {

                // add char
                b'A' | b'C' | b'G' | b'T' => {
                    let c: u8 = *c;

                    let mut out_edges = EMPTY_EDGE_LIST;
                    let mut in_edges = EMPTY_EDGE_LIST;

                    if (counter as usize) >= edges.len() {
                        out_edges = edges[counter as usize].0.clone();
                        in_edges = edges[counter as usize].1.clone();
                    }

                    match current_index_letter_type {
                        SolidString => {
                            if prev_index_letter_type == current_index_letter_type {
                                out_edges.push(counter as u32 +1);
                            }

                            // add prev edge
                            if counter > 0 && prev_index_letter_type == current_index_letter_type {
                                in_edges.push(counter as u32 - 1);
                            }


                        },

                        DegenerateLetter => {
                            // add next edge
                            if prev_index_letter_type == current_index_letter_type {
                                // out_edges.push(counter as u32 +1);
                            }

                            // add prev edge
                            if counter > 0 && prev_index_letter_type == current_index_letter_type {
                                // in_edges.push(counter as u32 - 1);
                            }

                        },


                        _ => { current_index_letter_type = DegenerateLetter }
                    };



                    // where seeds end
                    // add this to start of next solid string
                    if prev_index_letter_type == DegenerateLetter && current_index_letter_type == EndElastic {
                        match seed_indices_stop.as_ref() {
                            Some(indices) => {
                                if prev_index_letter_type == DegenerateLetter && current_index_letter_type == SolidString {
                                    in_edges = EMPTY_EDGE_LIST;
                                }

                                for idx in indices {
                                    in_edges.push(*idx);
                                }

                                if prev_index_letter_type == DegenerateLetter && current_index_letter_type == SolidString {
                                    seed_indices_stop = None;
                                }

                            },
                            _ => {}
                        }
                    }



                    if (counter as usize) >= edges.len() {
                        edges[counter as usize] = (in_edges.clone(), out_edges.clone());
                    } else {
                        edges.push((in_edges.clone(), out_edges.clone()));
                    }

                    data.push(c);
                    counter += 1;


                    prev_index_letter_type = current_index_letter_type;
                },

                // add multiple out edges
                b'{' => {
                    if counter != 0 {
                        length += 1;
                    }

                    prev_index_letter_type = StartElastic;
                    current_index_letter_type = DegenerateLetter;
                    // save the index of the last char in solid string
                    last_solid_string_index = Some(prev as u32);
                    seed_indices_start = Some(vec![counter as u32]);
                    seed_indices_stop = Some(vec![]);
                },

                // add multiple in edges
                b'}' => {
                    length += 1;

                    // add the last stop index
                    if prev_index_letter_type != SeedBoundary {
                        edges.push(EMPTY_EDGE_PIAR);
                        seed_indices_stop.as_mut().map(|v| v.push(counter as u32));
                    }
                    // save the index of the last char in last seed
                    // let x = (counter + 1) as u32;

                    // update seed out edges
                    // next_solid_string_index = Some(x);
                    match seed_indices_stop.as_ref() {
                      Some(indices) => {
                            eprintln!("Stop indices {:?}", indices);
                            for idx in indices {
                                edges[*idx as usize].1.push(counter as u32);
                            }
                        }
                        _ => panic!("")
                    }

                    // update seed in edges
                    match seed_indices_start.as_ref() {
                        Some(indices) => {
                            eprintln!("Seed start indices {:?} {:?}", indices, last_solid_string_index);

                            for idx in indices {
                                edges[*idx as usize].0.push(last_solid_string_index.unwrap());
                            }
                        },
                        None => {}
                    }

                    // Update prev solid string out edges
                    // with seed starts
                    match last_solid_string_index {
                        Some(idx) => {
                            let idx = idx as usize;
                            seed_indices_start.map(|indices| {
                                for i in indices {
                                    edges[idx].1.push(i);
                                }
                            });
                        },
                        None => {} // started at dt
                    };

                    last_solid_string_index = None;
                    seed_indices_start = None;

                    current_index_letter_type = EndElastic;
                    prev_index_letter_type = DegenerateLetter;
                },

                // Seed boundary
                b',' => {


                    // accumulate seed starts
                    match seed_indices_start.as_mut() {
                        Some(v) => v.push(counter as u32),
                        // empty seed
                        None => panic!(""),
                    };

                    // accumulate seed ends
                    match seed_indices_stop.as_mut() {
                        Some(v) => v.push(counter as u32 - 1),
                        // empty seed
                        None => panic!(""),
                    };

                    current_index_letter_type = DegenerateLetter;
                    prev_index_letter_type = SeedBoundary;
                },

                _ => panic!("Malformed EDS {}", *c as char)
            }
        }

        EDT {
            data,
            edges,
            size: counter as u32,
            length,
        }
    }


}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_str() {

        // starts with solid string
        // let ed_string = "AT{TCC,C}AA";
        let ed_string = "AT{TCC,AG,C}AA";
        // let ed_string = "AT{TCC,C,}AA";
        let edt = EDT::from_str(ed_string);
        eprintln!("{}", ed_string);
        eprintln!("{:?}", edt.edges);
        eprintln!("{}", edt.data.len());

        /*
        // starts with solid string
        let ed_string = "ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATATTGAT{T,C}CCTA\
                         {T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGTTG{A,G}CTC{T,G}\
                         TGAGTGAGCTTGCGAGATA";
        let edt = EDT::from_str(ed_string);
        eprintln!("{:?}", edt);
        eprintln!("{}", edt.data.len());

        // starts with degenerate letter
        let ed_string = "{TCGA,CTA,A}ATCGATGGG{T,C}AACTT{T,G}AG{G,T}CCGGTTTATAT\
                         TGAT{T,C}CCTA{T,G}{T,A}{A,T}A{T,A}GGGGGTCCTTTGCTTGCTGT\
                         TG{A,G}CTC{T,G}TGAGTGAGCTTGCGAGATA";
        // let result = 2 + 2;
        // assert_eq!(result, 4);
         */
    }
}

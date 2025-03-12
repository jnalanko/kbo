use std::{cmp::min, ops::Range};

use sbwt::SbwtIndexVariant;

// Pads with dollars from the left if there is not full k-mer
fn get_kmer_ending_at(query: &[char], end_pos: usize, k: usize) -> Vec<char> {
	let mut query_kmer = Vec::<char>::new();
	if end_pos >= k-1 {
		query_kmer.extend(&query[end_pos + 1 - k .. end_pos+1]);
	} else {
		let n_dollars = -(end_pos as isize - k as isize + 1);
		assert!(n_dollars > 0);
		for _ in 0..n_dollars {
			query_kmer.push('$');
		}
		query_kmer.extend(&query[0..end_pos+1]);
	}
	assert!(query_kmer.len() == k);
	query_kmer
}

fn longest_common_prefix(x: &[u8], y: &[u8]) -> usize {
	let mut len = 0_usize;
	for i in 0..min(x.len(), y.len()) {
		if x[i] == y[i] {
			len += 1;
		} else {
			break;
		}
	}
	len
}

/*
fn longest_common_suffix(x: &[char], y: &[char]) -> usize {
	let mut len = 0_usize;
	for i in 0..min(x.len(), y.len()) {
		if x[x.len() - 1 - i] == y[y.len() - 1 - i] {
			// The if-check will not go to negatives because of the constraint on i
			len += 1;
		} else {
			break;
		}
	}
	len
}
*/

fn chars_to_bytes(chars: Vec<char>) -> Vec<u8> {
    chars.iter().flat_map(|&c| c.to_string().into_bytes()).collect()
}

enum Variant {
	Subsitution(char), // Substitution to this character
	Deletion(usize), // How many
	Insertion(Vec<char>), // Insertion of these characters
}

struct VariantCalls {
	calls: Vec<(usize, Variant)> // Position, variant
}
#[allow(missing_docs)] // Will document when I know what this does
pub fn call_variants(
    sbwt: &SbwtIndexVariant,
    query: &[char],
    translation: &[char],
    original_ms: &[(usize, Range<usize>)],
    derand_ms: &[i64],
    k: usize,
) -> VariantCalls {

	let mut calls: Vec<(usize, Variant)> = vec![];

	let sbwt = match sbwt {
		SbwtIndexVariant::SubsetMatrix(index) => index,
		_ => panic!("Only SbwtIndexVariant::SubsetMatrix is supported"),
	};

	for (i, dms) in derand_ms.iter().enumerate() {
		let dms = derand_ms[i];
		if dms == 0 && translation[i] != 'M' {
			// Go to closest unique match position to the right
			let mut unique_pos: Option<(usize, usize)> = None; // (query pos, colex rank)
			for j in i+1..min(i+k+1, query.len()) {
				if original_ms[j].1.len() == 0 {
					unique_pos = Some((j, original_ms[j].1.start));
					break;
				}
			}

			if let Some((unique_end_pos, colex)) = unique_pos {
				assert!(unique_end_pos >= i);
				let unique_match_len = unique_end_pos - i;
				let ref_kmer = sbwt.access_kmer(colex);
				let query_kmer = chars_to_bytes(get_kmer_ending_at(query, unique_end_pos, k));

				if longest_common_prefix(&ref_kmer, &query_kmer) + unique_match_len == k-1 {
					// Single nucleotide substitution at i
					calls.push((i, Variant::Subsitution(query[i])));
				}
				else if longest_common_prefix(&ref_kmer, &query_kmer[1..]) + unique_match_len >= k-1 {
					// Single nucleotide deletion at i
					
					// We check for >= k-1 instead of == k-1 because we can have random matches after
					// the deletion point. For example:

					//           deleted character
					//                  v
					// Ref   ACACGCTAGCAGGCTGACTCGAT
					//        \\\\\\\\\\\|||||||||||
					// Query GACACGCTAGCAGCTGACTCGAT
					//
					// Here the deleted G happens to match the next character
					// in the query after the deletion point.

					calls.push((i, Variant::Deletion(1)));
				}

				else if longest_common_prefix(&ref_kmer[1..], &query_kmer) + unique_match_len >= k-1 {
					// Single nucleotide insertion at i
					// We check for >= k-1 instead of == k-1 for the same reason as in deletions

					calls.push((i, Variant::Insertion(1)));
				}
			}
		}
	}

    VariantCalls{calls}
}
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

#[derive(Debug)]
enum Variant {
	Subsitution(char), // Substitution to this character
	Deletion(usize), // How many
	Insertion(Vec<char>), // Insertion of these characters
}

#[derive(Debug)]
struct VariantCalls {
	calls: Vec<(usize, Variant)> // Position, variant
}

#[allow(missing_docs)] // Will document when I know what this does
fn call_variants(
    sbwt: &SbwtIndexVariant,
    query: &[char],
    original_ms: &[(usize, Range<usize>)],
	significant_match_threshold: usize,
    derand_ms: &[i64],
    k: usize,
) -> VariantCalls {

	let d = significant_match_threshold; // Shorthand

	let mut calls: Vec<(usize, Variant)> = vec![];

	let sbwt = match sbwt {
		SbwtIndexVariant::SubsetMatrix(index) => index,
		_ => panic!("Only SbwtIndexVariant::SubsetMatrix is supported"),
	};

    let mut prev_dms = 0_i64;
	for (i, &dms) in derand_ms.iter().enumerate() {
		if dms == 0 || (dms < k as i64 && prev_dms == k as i64) {
			// Go to closest unique match position to the right
			let mut unique_pos: Option<(usize, usize, usize)> = None; // (query pos, ms, colex rank)
			eprintln!("{} {:?} {} {}", i, unique_pos, i+k+1, query.len());
			for j in i+1..min(i+k+1, query.len()) {
				dbg!(&j, &original_ms[j]);
				if original_ms[j].1.len() == 1 {
					unique_pos = Some((j, original_ms[j].0, original_ms[j].1.start));
					break;
				}
			}

			if let Some((unique_end_pos, unique_match_len, colex)) = unique_pos {
				assert!(unique_end_pos >= i);
				let ref_kmer = sbwt.access_kmer(colex);
				let query_kmer = chars_to_bytes(get_kmer_ending_at(query, unique_end_pos, k));

				eprintln!("{}", String::from_utf8_lossy(&ref_kmer));
				eprintln!("{}", String::from_utf8_lossy(&query_kmer));
				eprintln!("{} {}", longest_common_prefix(&ref_kmer, &query_kmer[1..]), unique_match_len);

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

					calls.push((i, Variant::Insertion(vec![query[i]])));
				}
			}
		}
        prev_dms = dms;
	}

    VariantCalls{calls}
}


#[cfg(test)]
mod tests {

    use crate::{build, derandomize::derandomize_ms_vec, index::{query_sbwt, BuildOpts}, translate::translate_ms_vec};

    use super::*;

    #[test]
    fn test_variant_calling() {

		let k = 7;

        //                                 deleted character
        //                                        v
        let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATAC";
        let query =                b"GACACGCTAGCAGCTGACTCGAT";

        let (sbwt, lcs) = build(&[reference.to_vec()], BuildOpts{ build_select: true, k, ..Default::default() });
		let n_kmers = match &sbwt {
			SbwtIndexVariant::SubsetMatrix(s) => s.n_kmers(),
			_ => panic!("Only plain matrix sbwt is supported"),
		};

		let threshold = crate::derandomize::random_match_threshold(k, n_kmers, 4_usize, 0.001_f64);

        let noisy_ms = query_sbwt(query, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);
        let translated = translate_ms_vec(&derand_ms, k, threshold);

        eprintln!("{:?}", noisy_ms);
        eprintln!("{:?}", derand_ms);
        eprintln!("{:?}", translated);

        let query_chars: Vec<char> = query.to_vec().iter().map(|c| *c as char).collect();
        let variants = call_variants(&sbwt, &query_chars, &noisy_ms, threshold, &derand_ms, k);

        dbg!(variants);

        assert!(false);
    }

}
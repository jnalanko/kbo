use std::{cmp::min, ops::Range, process::exit};

use sbwt::{ContractLeft, ExtendRight, LcsArray, SbwtIndex, SbwtIndexVariant, StreamingIndex, SubsetMatrix};

// Pads with dollars from the left if there is not full k-mer
fn get_kmer_ending_at(query: &[u8], end_pos: usize, k: usize) -> Vec<u8> {
	let mut query_kmer = Vec::<u8>::new();
	if end_pos >= k-1 {
		query_kmer.extend(&query[end_pos + 1 - k .. end_pos+1]);
	} else {
		let n_dollars = -(end_pos as isize - k as isize + 1);
		assert!(n_dollars > 0);
		for _ in 0..n_dollars {
			query_kmer.push(b'$');
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

fn longest_common_suffix(x: &[u8], y: &[u8]) -> usize {
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

fn chars_to_bytes(chars: Vec<char>) -> Vec<u8> {
    chars.iter().flat_map(|&c| c.to_string().into_bytes()).collect()
}

#[derive(Debug, Eq, PartialEq)]
struct Variant {
	query_pos: usize,
	query_chars: Vec<u8>, // If empty, it's a query deletion
	ref_chars: Vec<u8>, // If empty, it's a query insertion 
}

// Assumes things about the inputs: todo: document
fn get_variant_length(ms: &[(usize, Range<usize>)], common_suffix_len: usize, significant_match_threshold: usize) -> usize {
	assert!(ms.len() > common_suffix_len);

	// Todo: make sure that none of the arithmetic can result in negative values

	let mismatch_pos = ms.len() - common_suffix_len - 1;
	// Go to the right until we are that match threshold again
	for i in mismatch_pos..ms.len() {
		let match_len = ms[i].0;
		//dbg!(i,mismatch_pos,match_len,common_suffix_len);
		if match_len >= significant_match_threshold {
			return i - mismatch_pos + 1 - match_len;
		}
	}

	// There is no significant match to the right -> we use the
	// rightmost match.
	let match_len = ms.last().unwrap().0;
	//dbg!(match_len, ms.len(), mismatch_pos);
	(ms.len() - 1) - mismatch_pos + 1 - match_len
}

fn get_rightmost_significant_peak(ms: &[(usize, Range<usize>)], significant_match_threshold: usize) -> Option<usize> {
	assert!(!ms.is_empty());
	for i in (0..ms.len()-1).rev() {
		let here = ms[i].0;
		let next = ms[i+1].0;
		if here >= significant_match_threshold && here > next {
			return Some(i);
		}
	}
	None
}

#[allow(clippy::too_many_arguments)]
fn resolve_variant(
	query_kmer: &[u8], 
	ref_kmer: &[u8], 
	ms_vs_query: &[(usize, Range<usize>)], // Slice of length k
	ms_vs_ref: &[(usize, Range<usize>)], // Slice of length k
	common_suffix_len: usize,
	k: usize,
	significant_match_threshold: usize,
) -> Option<(Vec<u8>, Vec<u8>)> { // Returns (query chars, ref chars)

	assert!(ms_vs_query.len() == k);
	assert!(ms_vs_ref.len() == k);
	assert!(common_suffix_len > 0);
	//assert!(unique_end_pos >= i);

	eprintln!("{:?}", ms_vs_ref);
	eprintln!("{:?}", ms_vs_query);

	let query_ms_peak = get_rightmost_significant_peak(ms_vs_ref, significant_match_threshold);
	let ref_ms_peak = get_rightmost_significant_peak(ms_vs_query, significant_match_threshold);

	dbg!(query_ms_peak, ref_ms_peak);

	if let (Some(query_ms_peak), Some(ref_ms_peak)) = (query_ms_peak, ref_ms_peak) {
		let suffix_match_start = k - common_suffix_len;

		// Negative gap means overlap 
		let query_gap = suffix_match_start as isize - query_ms_peak as isize - 1;
		let ref_gap = suffix_match_start as isize - ref_ms_peak as isize - 1;

		dbg!(query_gap, ref_gap);

		if query_gap > 0 && ref_gap > 0 {
			return Some(
				(
					query_kmer[query_ms_peak+1..suffix_match_start].to_vec(), // Query chars
					ref_kmer[ref_ms_peak+1..suffix_match_start].to_vec() // Ref chars
				)
			)
		} else {
			let query_overlap = -query_gap;
			let ref_overlap = -ref_gap;
			assert!(query_overlap != ref_overlap); // This is impossible?
			let variant_len = (query_overlap - ref_overlap).unsigned_abs();
			if query_overlap > ref_overlap {
				// Deletion in query
				return Some(
					(
						vec![], // Query chars
						ref_kmer[ref_ms_peak + 1 .. ref_ms_peak + 1 + variant_len].to_vec() // Ref chars
					)
				);
			}  else {
				// Insertion in query
				return Some(
					(
						query_kmer[query_ms_peak + 1 .. query_ms_peak + 1 + variant_len].to_vec(), // Query chars
						vec![] // Ref chars
					)
				);
			}
		}

	}
	
	None // Could not resolve variant
}

#[allow(missing_docs)] // Will document when I know what this does
fn call_variants(
	sbwt_ref: &SbwtIndex<SubsetMatrix>,
	lcs_ref: &LcsArray,
	sbwt_query: &SbwtIndex<SubsetMatrix>,
	lcs_query: &LcsArray,
    query: &[u8],
	significant_match_threshold: usize,
    k: usize,
) -> Vec<Variant> {

	let d = significant_match_threshold; // Shorthand

	let mut calls: Vec<Variant> = vec![];

	let index_ref = StreamingIndex::new(&sbwt_ref, &lcs_ref);
	let index_query = StreamingIndex::new(&sbwt_query, &lcs_query);
	let ms_vs_ref = index_ref.matching_statistics(&query);

	eprintln!("{:?}", ms_vs_ref.iter().map(|x| x.0).collect::<Vec::<usize>>());

	for i in 1..query.len() {
		if ms_vs_ref[i].0 < ms_vs_ref[i-1].0 && ms_vs_ref[i-1].0 >= d && ms_vs_ref[i].0 < d {
			// Go to closest unique match position to the right
			eprintln!("{} {} {}", i, i+k+1, query.len());
			for j in i+1..min(i+k+1, query.len()) {
				if ms_vs_ref[j].0 >= significant_match_threshold && ms_vs_ref[j].1.len() == 1 {
					eprintln!("Investigating positions {} {}, ms[j] = {}", i, j, ms_vs_ref[j].0);
					let ref_colex = ms_vs_ref[j].1.start;

					let query_kmer = get_kmer_ending_at(query, j, k);
					let ref_kmer = sbwt_ref.access_kmer(ref_colex);
					let suffix_match_len = longest_common_suffix(&query_kmer, &ref_kmer);

					eprintln!("{}", String::from_utf8_lossy(&ref_kmer));
					eprintln!("{}", String::from_utf8_lossy(&query_kmer));
					dbg!(suffix_match_len);
					
					// MS vectors for k-mers (todo: use slices of global MS vector?)
					let ms_vs_ref = index_ref.matching_statistics(&query_kmer);
					let ms_vs_query = index_query.matching_statistics(&ref_kmer);

					if let Some(var) = resolve_variant(&query_kmer, &ref_kmer, &ms_vs_query, &ms_vs_ref, suffix_match_len, k, significant_match_threshold) {
						calls.push(Variant{query_chars: var.0, ref_chars: var.1, query_pos: i});
					}
					break;
				}
			}


		}
	}

	calls
}


#[cfg(test)]
mod tests {

    use std::cmp::max;

    use sbwt::{BitPackedKmerSortingMem, SbwtConstructionAlgorithm, SbwtIndexBuilder};

    use crate::{build, derandomize::derandomize_ms_vec, index::{query_sbwt, BuildOpts}, translate::translate_ms_vec};

    use super::*;

	fn run_variant_calling(query: &[u8], reference: &[u8], k: usize) -> Vec<Variant> {
        let (sbwt_ref, lcs_ref) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[reference]);
        let (sbwt_query, lcs_query) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[query]);

		let threshold = crate::derandomize::random_match_threshold(k, max(sbwt_ref.n_kmers(), sbwt_query.n_kmers()), 4_usize, 0.001_f64);

        eprintln!("t = {}", threshold);

		call_variants(&sbwt_ref, lcs_ref.as_ref().unwrap(), &sbwt_query, lcs_query.as_ref().unwrap(), query, threshold, k)

	}

	#[test]
	fn variant_call_test_single_base_substitution() {
		let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
		let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCTTCAACGTG";
		//                                                                 ^ here

		let variants = run_variant_calling(query, reference, 20);
		dbg!(&variants);

		assert_eq!(variants, vec![Variant{query_pos: 49, ref_chars: vec![b'A'], query_chars: vec![b'T']}]);
	}

	#[test]
	fn variant_call_test_multi_base_substitution() {
		//                                             **
		let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
		let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAGCGTCTATTGTACCAATCGGCATCAACGTG";
		//                                             ***

		let variants = run_variant_calling(query, reference, 30);
		dbg!(&variants);

		assert_eq!(variants, vec![Variant{query_pos: 29, ref_chars: b"AA".to_vec(), query_chars: b"GCG".to_vec()}]);
	}

	#[test]
	fn variant_call_test_single_base_insertion_case1() {
		// Case: Non-overlapping reference intervals

		let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
		let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCAGTCAACGTG";
		//                                                                  ^ here

		let variants = run_variant_calling(query, reference, 20);
		dbg!(&variants);

		assert_eq!(variants, vec![Variant{query_pos: 50, ref_chars: vec![], query_chars: vec![b'G']}]);
	}

	#[test]
	fn variant_call_test_single_base_insertion_case2() {
		// Case: Overlapping reference intervals

		let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
		let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCAATCAACGTG";
		//                                                                  ^ here

		let variants = run_variant_calling(query, reference, 20);
		dbg!(&variants);

		assert_eq!(variants, vec![Variant{query_pos: 50, ref_chars: vec![], query_chars: vec![b'A']}]);
	}

	#[test]
	fn variant_call_test_single_base_deletion_case1() {
		// Case: Non-overlapping query intervals

		let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCAGTCAACGTG";
		let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
		//                                                                  ^ here

		let variants = run_variant_calling(query, reference, 20);
		dbg!(&variants);

		assert_eq!(variants, vec![Variant{query_pos: 50, ref_chars: vec![b'G'], query_chars: vec![]}]);
	}

	#[test]
	fn variant_call_test_single_base_deletion_case2() {
		// Case: overlapping query intervals

		let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATTCAACGTG";
		let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
		//                                                                  ^ here

		let variants = run_variant_calling(query, reference, 20);
		dbg!(&variants);

		assert_eq!(variants, vec![Variant{query_pos: 51, ref_chars: vec![b'T'], query_chars: vec![]}]);
	}


    #[test]
    fn variant_call_test_variants_in_same_query() {
        //                                 deleted character    substituted        inserted
        //                                        v                 v                v
        let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATACTATGTGTTATAGCAATTCGGATCGATCGA";
        let query =      b"TCGTGGATCGATACACGCTAGCAGCTGACTCGATGGGATACCATGTGTTATAGCAATTCCGGATCGATCGA";

		let variants = run_variant_calling(query, reference, 20);
		dbg!(&variants);

		assert_eq!(variants[0], Variant{query_pos: 24, query_chars: vec![], ref_chars: vec![b'G']});
		assert_eq!(variants[1], Variant{query_pos: 41, query_chars: vec![b'C'], ref_chars: vec![b'T']});
		assert_eq!(variants[2], Variant{query_pos: 59, query_chars: vec![b'C'], ref_chars: vec![]});
		assert_eq!(variants.len(), 3);
    }

}
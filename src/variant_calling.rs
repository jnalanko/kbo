use std::{cmp::min, ops::Range};

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
	Substitution((u8, u8)), // Substitution from this to that character
	Deletion(Vec<u8>), // Deletion of these characters 
	Insertion(Vec<u8>), // Insertion of these characters
}

#[derive(Debug)]
struct VariantCalls {
	calls: Vec<(usize, Variant)> // Position, variant
}

// Assumes things about the inputs: todo: document
fn get_variant_length(ms: &[(usize, Range<usize>)], common_suffix_len: usize, significant_match_threshold: usize) -> usize {
	assert!(ms.len() > common_suffix_len);

	// Todo: make sure that none of the arithmetic can result in negative values

	let mismatch_pos = ms.len() - common_suffix_len - 1;
	// Go to the right until we are that match threshold again
	for i in mismatch_pos..ms.len() {
		let match_len = ms[i].0;
		if match_len >= significant_match_threshold {
			dbg!(i, mismatch_pos, match_len, common_suffix_len, ms.len());
			return i - mismatch_pos + 1 - match_len;
		}
	}

	// There is no significant match to the right -> we use the
	// rightmost match.
	let match_len = ms.last().unwrap().0;
	(ms.len() - 1) - mismatch_pos + 1 - match_len
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
) -> Option<Variant> {

	assert!(ms_vs_query.len() == k);
	assert!(ms_vs_ref.len() == k);
	assert!(common_suffix_len > 0);
	//assert!(unique_end_pos >= i);

	let variant_end = k - common_suffix_len; // Exclusive end

	let query_var_len = get_variant_length(ms_vs_ref, common_suffix_len, significant_match_threshold);
	let ref_var_len = get_variant_length(ms_vs_query, common_suffix_len, significant_match_threshold);

	assert!(variant_end > query_var_len);
	assert!(variant_end > ref_var_len);

	if query_var_len == 0 && ref_var_len > 0 {
		// Deletion in query
		let seq = ref_kmer[variant_end - ref_var_len .. variant_end].to_vec();
		return Some(Variant::Deletion(seq));
	} else if query_var_len > 0 && ref_var_len == 0 {
		// Insertion in query
		let seq = query_kmer[variant_end - query_var_len .. variant_end].to_vec(); // Todo: figure out indices
		return Some(Variant::Insertion(seq));
	} else if query_var_len == 1 && ref_var_len == 1 {
		// Substitution
		let ref_char = ref_kmer[variant_end-1];
		let query_char = query_kmer[variant_end-1];
		return Some(Variant::Substitution((ref_char, query_char)));
	}

	// Todo: Check bases before the variation site
	
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
) -> VariantCalls {

	let d = significant_match_threshold; // Shorthand

	let mut calls: Vec<(usize, Variant)> = vec![];

	let index_ref = StreamingIndex::new(&sbwt_ref, &lcs_ref);
	let index_query = StreamingIndex::new(&sbwt_query, &lcs_query);
	let ms_vs_ref = index_ref.matching_statistics(&query);

	eprintln!("{:?}", ms_vs_ref.iter().map(|x| x.0).collect::<Vec::<usize>>());

	for i in 1..query.len() {
		if ms_vs_ref[i].0 < ms_vs_ref[i-1].0 && ms_vs_ref[i-1].0 >= d && ms_vs_ref[i].0 < d {
			// Go to closest unique match position to the right
			eprintln!("{} {} {}", i, i+k+1, query.len());
			for j in i+1..min(i+k+1, query.len()) {
				if ms_vs_ref[j].1.len() == 1 {
					eprintln!("Investigating positions {} {}", i, j);
					let ref_colex = ms_vs_ref[j].1.start;
					let common_suffix_len = min(ms_vs_ref[j].0, j-i); // Don't go over the variant position

					let query_kmer = get_kmer_ending_at(query, j, k);
					let ref_kmer = sbwt_ref.access_kmer(ref_colex);

					eprintln!("{}", String::from_utf8_lossy(&ref_kmer));
					eprintln!("{}", String::from_utf8_lossy(&query_kmer));

					// MS vectors for k-mers (todo: use slices of global MS vector?)
					let ms_vs_ref = index_ref.matching_statistics(&query_kmer);
					let ms_vs_query = index_query.matching_statistics(&ref_kmer);

					if let Some(var) = resolve_variant(&query_kmer, &ref_kmer, &ms_vs_query, &ms_vs_ref, common_suffix_len, k, significant_match_threshold) {
						calls.push((i, var));
					}
					
					dbg!(&calls);
				}
			}


		}
	}

    VariantCalls{calls}
}


#[cfg(test)]
mod tests {

    use std::cmp::max;

    use sbwt::{BitPackedKmerSortingMem, SbwtConstructionAlgorithm, SbwtIndexBuilder};

    use crate::{build, derandomize::derandomize_ms_vec, index::{query_sbwt, BuildOpts}, translate::translate_ms_vec};

    use super::*;

    #[test]
    fn test_variant_calling() {

		let k = 20;

        //                                 deleted character    substituted        inserted
        //                                        v                 v                v
        let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATACTATGTGTTATAGCAATTCGGATCGATCGA";
        let query =      b"TCGTGGATCGATACACGCTAGCAGCTGACTCGATGGGATACCATGTGTTATAGCAATTCCGGATCGATCGA";


        let (sbwt_ref, lcs_ref) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[reference]);
        let (sbwt_query, lcs_query) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[query]);

		let threshold = crate::derandomize::random_match_threshold(k, max(sbwt_ref.n_kmers(), sbwt_query.n_kmers()), 4_usize, 0.001_f64);

        eprintln!("t = {}", threshold);

		let variants = call_variants(&sbwt_ref, lcs_ref.as_ref().unwrap(), &sbwt_query, lcs_query.as_ref().unwrap(), query, threshold, k);

        dbg!(variants);

        assert!(false);
    }

}
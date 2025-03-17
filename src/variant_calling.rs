use std::{cmp::min, ops::Range};

use sbwt::{ContractLeft, ExtendRight, SbwtIndexVariant, StreamingIndex};

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
	Subsitution((u8, u8)), // Substitution from this to that character
	Deletion(Vec<u8>), // Deletion of these characters 
	Insertion(Vec<u8>), // Insertion of these characters
}

#[derive(Debug)]
struct VariantCalls {
	calls: Vec<(usize, Variant)> // Position, variant
}

fn locate(colex: usize) -> usize {
	todo!();
}

fn get_variant_length(ms: &[(usize, Range<usize>)], significant_match_threshold: usize) -> usize {
	// Go to the right until we are that match threshold again
	todo!();
}

#[allow(clippy::too_many_arguments)]
fn resolve_variant(
	query_kmer_end: usize, // Text position
	ref_kmer_end: usize, // Text position
	query: &[u8], 
	reference: &[u8], 
	ms_vs_query: &[(usize, Range<usize>)],
	ms_vs_ref: &[(usize, Range<usize>)],
	k: usize,
	significant_match_threshold: usize,
) -> Option<Variant> {

	//assert!(unique_end_pos >= i);
	let query_kmer = get_kmer_ending_at(&query, query_kmer_end, k);
	let ref_kmer = get_kmer_ending_at(&reference, ref_kmer_end, k);

	eprintln!("{}", String::from_utf8_lossy(&ref_kmer));
	eprintln!("{}", String::from_utf8_lossy(&query_kmer));

	let query_var_length = get_variant_length(ms_vs_ref, significant_match_threshold);
	let ref_var_length = get_variant_length(ms_vs_query, significant_match_threshold);

	if query_var_length == 0 && ref_var_length > 0 {
		// Deletion in query
		let seq = ref_kmer[..].to_vec(); // Todo: figure out indices
		return Some(Variant::Deletion(seq));
	} else if (query_var_length > 0 && ref_var_length == 0) {
		// Insertion in query
		let seq = query_kmer[..].to_vec(); // Todo: figure out indices
		return Some(Variant::Insertion(seq));
	} else if (query_var_length == 1 && ref_var_length == 1) {
		// Substitution
		let ref_char = ref_kmer[0]; // Todo: figure out indices
		let query_char = query_kmer[0]; // Todo: figure out indices
		return Some(Variant::Subsitution((ref_char, query_char)));
	}

	// Todo: Check bases before the variation site
	
	return None; // Could not resolve variant
}

#[allow(missing_docs)] // Will document when I know what this does
fn call_variants<E: ExtendRight, C: ContractLeft>(
    index_ref: &StreamingIndex<E, C>,
    index_query: &StreamingIndex<E, C>,
	reference: &[u8],
    query: &[u8],
	significant_match_threshold: usize,
    k: usize,
) -> VariantCalls {

	let d = significant_match_threshold; // Shorthand

	let mut calls: Vec<(usize, Variant)> = vec![];

	let ms_vs_ref = index_ref.matching_statistics(&query);
	let ms_vs_query = index_ref.matching_statistics(&reference);

    let mut prev_dms = 0_i64;
	for i in 1..query.len() {
		if ms_vs_ref[i].0 < ms_vs_ref[i-1].0 && ms_vs_ref[i-1].0 >= d {
			// Go to closest unique match position to the right
			let mut unique_pos: Option<(usize, usize, usize)> = None; // (query pos, ms, colex rank)
			eprintln!("{} {:?} {} {}", i, unique_pos, i+k+1, query.len());
			for j in i+1..min(i+k+1, query.len()) {
				if ms_vs_ref[j].1.len() == 1 {
					let query_pos = j;
					let ref_colex = ms_vs_ref[j].1.start;
					let ref_pos = locate(ref_colex);
					resolve_variant(j, ref_pos, &query, &reference, &ms_vs_query, &ms_vs_ref, k, d);
				}
			}


		}
	}

    VariantCalls{calls}
}


#[cfg(test)]
mod tests {

    use crate::{build, derandomize::derandomize_ms_vec, index::{query_sbwt, BuildOpts}, translate::translate_ms_vec};

    use super::*;

	/*
    #[test]
    fn test_variant_calling() {

		let k = 20;

        //                                 deleted character    substituted        inserted
        //                                        v                 v                v
        let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATACTATGTGTTATAGCAATTCGGATCGATCGA";
        let query =      b"TCGTGGATCGATACACGCTAGCAGCTGACTCGATGGGATACCATGTGTTATAGCAATTCCGGATCGATCGA";


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
        eprintln!("t = {}", threshold);

        let query_chars: Vec<char> = query.to_vec().iter().map(|c| *c as char).collect();
        let variants = call_variants(&sbwt, &query_chars, &noisy_ms, threshold, &derand_ms, k);

        dbg!(variants);

        assert!(false);
    }
	*/

}
// kbo: Spectral Burrows-Wheeler transform accelerated local alignment search
//
// Copyright 2024 Tommi Mäklin [tommi@maklin.fi].

// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.

// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
//! Derandomizing noisy _k_-bounded matching statistics.

use embed_doc_image::embed_doc_image;

/// Controls parameters of the derandomization algorithm.
///
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct DerandomizeOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
}

impl Default for DerandomizeOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::derandomize::DerandomizeOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// # let expected = kbo::derandomize::DerandomizeOpts::default();
    /// # assert_eq!(opts.max_error_prob, expected.max_error_prob);
    /// ```
    ///
    fn default() -> DerandomizeOpts {
        DerandomizeOpts {
            max_error_prob: 0.0000001,
        }
    }
}

/// Evaluates the CDF of _k_-bounded matching statistics random match distribution.
///
/// Computes the log-probability that a matching statistic with value
/// `t` or less that is the result of mapping a _k_-mer with
/// `alphabet_size` possible characters against an index containing
/// `n_kmers` _k_-mers was generated by chance.
///
/// This probability is given by the cumulative distribution function (CDF):
///
/// ![Formula: (1 - q^{t + 1})^n][latex_rmp],
///
/// where ![q = 1/s][latex_1os] with _s_ possible characters
/// `alphabet_size` and _n_ _k_-mers in the index
/// `n_kmers`. Derivation is given below.
///
/// # Distribution of prefix match lengths in random strings
///
/// Let _X_ be a random variable denoting the length _t_ of the
/// longest common prefix of two uniformly random and infinitely long
/// strings. If &nbsp; ![0 < p < 1][latex_p_bounds] &nbsp; is the probability that any pair of two
/// characters in the two strings are mismatched, then _X_ follows the
/// [geometric
/// distribution](https://en.wikipedia.org/wiki/Geometric_distribution)
/// with cumulative distribution function
///
/// ![P(X less or equal to t) = 1 - (1 - p)^{t + 1} with t a natural number][latex_geom_dist_cdf]
///
/// Suppose that an index contains _n_ uniformly random strings of
/// length _k_. Assuming that _k_ is large enough that the probability
/// that two strings from the index match by chance is neglible, the
/// longest common prefix of these two strings is reasonably
/// approximated by the distribution of the random variable _X_.
///
/// Now, let _M_ be a random variable denoting the length of
/// the longest common prefix between some string of length
/// _k_ and the entire index. Since the longest common
/// prefix of two pairs of strings is given by _X_,
/// _M_ is the maximum of _n_ independent random
/// variables with the same distribution as _X_ given by ![M
/// = max{X_1,...,X_n}][latex_M_max]. Because the variables
/// _X_ were assumed independent, the cumulative
/// distribution function of _M_ is their product
///
/// ![P(M less or equal to t) = P(X less or equal to t)^n = (1 - (1 - p)^{t + 1})^n][latex_M_cdf_1]
///
/// By noting that ![p = 1 - q][latex_p_to_q], we get the cumulative distribution function
///
/// ![P(M less or equal to t) = (1 - q^{t + 1})^n][latex_M_cdf_2]
///
/// Credit to [Jarno N. Alanko](https://jnalanko.net/) for deriving the random match distribution.
///
/// # Examples
/// ```rust
/// # use assert_approx_eq::assert_approx_eq;
/// use kbo::derandomize::log_rm_max_cdf;
///
/// let alphabet_size = 4;
/// let n_kmers = 20240921;
///
/// let res = log_rm_max_cdf(10, alphabet_size, n_kmers);
/// // `res` is -4.825812199808644
/// # assert_approx_eq!(res, -4.825812199808644, 1e-8);
/// ```
///
#[embed_doc_image("latex_rmp", "doc/img/latex_random_match_prob.svg")]
#[embed_doc_image("latex_p_bounds", "doc/img/latex_p_bounds.svg")]
#[embed_doc_image("latex_q_bounds", "doc/img/latex_q_bounds.svg")]
#[embed_doc_image("latex_geom_dist_cdf", "doc/img/latex_geom_dist_cdf.svg")]
#[embed_doc_image("latex_1os", "doc/img/latex_1_over_s.svg")]
#[embed_doc_image("latex_M_max", "doc/img/latex_M_max.svg")]
#[embed_doc_image("latex_M_cdf_1", "doc/img/latex_M_cdf_1.svg")]
#[embed_doc_image("latex_M_cdf_2", "doc/img/latex_M_cdf_2.svg")]
#[embed_doc_image("latex_p_to_q", "doc/img/latex_p_to_q.svg")]
pub fn log_rm_max_cdf(
    t: usize,
    alphabet_size: usize,
    n_kmers: usize,
) -> f64 {
    assert!(n_kmers > 0);
    assert!(alphabet_size > 0);

    n_kmers as f64 * (- ((1.0_f64.ln() - (alphabet_size as f64).ln()).exp()).powi(t as i32 + 1)).ln_1p()
}

/// Determines a lower bound for non-random _k_-bounded matching statistic values.
///
/// Computes the probabilities that the possible values for the
/// _k_-bounded matching statistics (MS) of a _k_-mer with size `k`
/// mapped against an index with `n_kmers` total _k_-mers and
/// `alphabet_size` possible values at each character are random
/// matches. Computation terminates when the MS value that produces a
/// random match probability below `max_error_prob` is found and
/// returned.
///
/// If no MS value passes the check, the function returns `k` instead.
///
/// # Examples
/// ```rust
/// use kbo::derandomize::random_match_threshold;
///
/// let k = 31;
/// let n_kmers = 20240921;
/// let alphabet_size = 4;
/// let max_error_prob = 0.01_f64;
///
/// let threshold = random_match_threshold(k, n_kmers, alphabet_size, max_error_prob);
/// // `threshold` is 15
/// # assert_eq!(threshold, 15);
/// ```
pub fn random_match_threshold(
    k: usize,
    n_kmers: usize,
    alphabet_size: usize,
    max_error_prob: f64,
) -> usize {
    assert!(k > 0);
    assert!(n_kmers > 0);
    assert!(alphabet_size > 0);
    assert!(max_error_prob <= 1_f64);
    assert!(max_error_prob > 0_f64);

    for i in 1..k {
	if log_rm_max_cdf(i, alphabet_size, n_kmers) > (-max_error_prob).ln_1p() {
	    return i;
	}
    }
    k
}

/// Derandomizes a single noisy _k_-bounded matching statistic.
///
/// Derandomizes the `current_noisy_ms` matching statistic (MS) based
/// on the `next_derand_ms` value obtained from the output of this
/// function for the next noisy MS when read left-to-right, the
/// _k_-mer size `k`, and the `threshold` which specifies a lower
/// bound to consider the MS a non-random match.
///
/// Positive values of the output i64 value mean that i64 characters
/// from the beginning of the k-mer match the reference, ie. same as
/// the MS, while negative values denote distance from the last
/// character in the last _k_-mer that produced a match.
///
/// # Examples
/// ## Noisy MS has only matches
/// ```rust
/// use kbo::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Noisy MS         : 1,2,3,3,3
/// // Derandomized MS  : 1,2,3,3,3
/// // Testing this pos :     |
///
/// let derand_ms = derandomize_ms_val(3, 3, 2, 3);
/// // `derand_ms` is 3
/// # assert_eq!(derand_ms, 3);
/// ```
///
/// ## Noisy MS has only noise
/// ```rust
/// use kbo::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Noisy MS         :  0, 0, 2, 1,0
/// // Derandomized MS  : -4,-3,-2,-1,0
/// // Testing this pos :        |
///
/// let derand_ms = derandomize_ms_val(2, -1, 2, 3);
/// // `derand_ms` is -2
/// # assert_eq!(derand_ms, -2);
/// ```
///
/// ## Noisy MS is at beginning of a full _k_-mer match
/// ```rust
/// use kbo::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Noisy MS         : 1,2,3, 1,2
/// // Derandomized MS  : 1,2,3,-1,0
/// // Testing this pos :     |
///
/// let derand_ms = derandomize_ms_val(3, -1, 2, 3);
/// // `derand_ms` is 3
/// # assert_eq!(derand_ms, 3);
/// ```
///
/// ## Noisy MS is at beginning of a partial _k_-mer match
/// ```rust
/// use kbo::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 4, threshold = 2
/// //
/// // Noisy MS         : 1,2,3,-1,0,1,2,3,4,4
/// // Derandomized MS  : 1,2,3,-1,0,1,2,3,4,4
/// // Testing this pos :     |
///
/// let derand_ms = derandomize_ms_val(3, -1, 2, 4);
/// // `derand_ms` is 3
/// # assert_eq!(derand_ms, 3);
/// ```
///
pub fn derandomize_ms_val(
    curr_noisy_ms: usize,
    next_derand_ms: i64,
    threshold: usize,
    k: usize,
) -> i64 {
    assert!(k > 0);
    assert!(threshold > 1);
    assert!(curr_noisy_ms <= k);
    assert!(next_derand_ms <= k as i64);

    // Default is to decrease MS by 1.
    let mut run: i64 = next_derand_ms - 1;

    if curr_noisy_ms == k {
	// Beginning of a full k-mer match
	run = k as i64;
    }

    if curr_noisy_ms > threshold && next_derand_ms < curr_noisy_ms as i64 {
	// Beginning of a partial k-mer match
	// Only useful if threshold > 1 and k > 3
	run = curr_noisy_ms as i64;
    }

    run
}

/// Derandomizes a sequence of noisy _k_-bounded matching statistics.
///
/// Iterates over a sequence of noisy _k_-bounded matching statistics
/// `ms` in reverse to identify values that are the result of random
/// matching between _k_-mers of size `k` and an index that the lower
/// bound `threshold` was calculated for.
///
/// # Examples
/// ```rust
/// use kbo::derandomize::derandomize_ms_vec;
///
/// let k = 3;
/// let threshold = 2;
/// let noisy_ms = vec![1,2,2,3,2,2,3,2,1,2,3,1,1,1,2,3,1,2];
///
/// let derand_ms = derandomize_ms_vec(&noisy_ms, k, threshold);
/// // `derand_ms` has [0,1,2,3,1,2,3,0,1,2,3,-1,0,1,2,3,-1,0]
/// # assert_eq!(derand_ms, vec![0,1,2,3,1,2,3,0,1,2,3,-1,0,1,2,3,-1,0]);
/// ```
///
pub fn derandomize_ms_vec(
    noisy_ms: &[usize],
    k: usize,
    threshold: usize,
) -> Vec<i64> {
    assert!(k > 0);
    assert!(threshold > 1);
    assert!(noisy_ms.len() > 2);

    let len = noisy_ms.len();
    let mut derand_ms: Vec<i64> = vec![0; len];

    // Traverse the matching statistics in reverse.
    derand_ms[len - 1] = if noisy_ms[len - 1] > threshold { noisy_ms[len - 1]} else { 0 } as i64;
    for i in 2..(len + 1) {
	derand_ms[len - i] = derandomize_ms_val(noisy_ms[len - i], derand_ms[len - i + 1], threshold, k);
    }

    derand_ms
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn log_rm_max_cdf() {
	let expected = vec![-1306319.1078024083,-318761.2492719044,-79220.9269610741,-19776.1823255263,-4942.2344281681,-1235.4454790664,-308.8543003470,-77.2131332649,-19.3032557026,-4.8258121998,-1.2064529421,-0.3016132288,-0.0754033068,-0.0188508267,-0.0047127067,-0.0011781767,-0.0002945442,-0.0000736360,-0.0000184090,-0.0000046023,-0.0000011506,-0.0000002876,-0.0000000719,-0.0000000180,-0.0000000045,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000];
	let alphabet_size = 4;
	let n_kmers = 20240921;
	let k = 1..32;
	k.for_each(|t| assert_approx_eq!(super::log_rm_max_cdf(t, alphabet_size, n_kmers), expected[t - 1], 1e-8f64));
    }

    #[test]
    fn random_match_threshold() {
	let expected = vec![15,18,22,25,28];
	let alphabet_size = 4;
	let n_kmers = 20240921;
	let k = 31;
	let factor = 1..6;
	factor.for_each(|i| assert_eq!(super::random_match_threshold(k, n_kmers, alphabet_size, (0.01_f64).powf(i as f64)), expected[i - 1]));
    }

    #[test]
    fn derandomize_ms_val_full_match() {
	// Parameters       : k = 3, threshold = 2
	//
	// Noisy MS         : 1,2,3,3,3
	// Derandomized MS  : 1,2,3,3,3
	// Testing this pos :     |

	let expected = 3;
	let got = super::derandomize_ms_val(3, 3, 2, 3);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_val_only_noise() {
	// Parameters       : k = 3, threshold = 2
	//
	// Noisy MS         :  0, 0, 2, 1,0
	// Derandomized MS  : -4,-3,-2,-1,0
	// Testing this pos :        |

	let expected = -2;
	let got = super::derandomize_ms_val(2, -1, 2, 3);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_val_beginning_of_full_match() {
	// Parameters       : k = 3, threshold = 2
	//
	// Noisy MS         : 1,2,3, 1,2
	// Derandomized MS  : 1,2,3,-1,0
	// Testing this pos :     |

	let expected = 3;
	let got = super::derandomize_ms_val(3, -1, 2, 3);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_val_beginning_of_partial_match() {
	// Parameters       : k = 4, threshold = 2
	//
	// Noisy MS         : 1,2,3,-1,0,1,2,3,4,4
	// Derandomized MS  : 1,2,3,-1,0,1,2,3,4,4
	// Testing this pos :     |

	let expected = 3;
	let got = super::derandomize_ms_val(3, -1, 2, 4);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_vec() {
	let noisy_ms = vec![1,2,2,3,2,2,3,2,1,2,3,1,1,1,2,3,1,2];
	let expected = vec![0,1,2,3,1,2,3,0,1,2,3,-1,0,1,2,3,-1,0];
	let got = super::derandomize_ms_vec(&noisy_ms, 3, 2);

	assert_eq!(got, expected);
    }
}

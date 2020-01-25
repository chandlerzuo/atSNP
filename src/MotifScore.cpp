#include "MotifScore.h"

/*
Compute likelihood ratio scores for SNPs' effect on motif matching.
@arg _motif_library The list object containing a 'matrix' component, which is a list of position weight matrices.
@arg _snpinfo A list object containing two components:
sequence_matrix: is a matrix for the sequences around each SNP. Each column corresponds to a SNP, and each row corresponds to a nucleotide position.
a1: A vector for the nucleobases on the reference genome at the SNP locations.
a2: A vector for the nucleobases after SNP.
@return A list of objects.
*/
extern "C" SEXP motif_score(SEXP _motif_library, SEXP _snpinfo) {
	//parse the _motif_library
	Rcpp::List pwms(_motif_library);

	//parse _snpinfo
	Rcpp::List snpinfo(_snpinfo);
	SEXP _sequence_matrix(snpinfo["sequence_matrix"]);
	IntegerMatrix sequence_matrix(_sequence_matrix);
	SEXP _ref_base_codes(snpinfo["ref_base"]);
	IntegerVector ref_base_codes(_ref_base_codes);
	SEXP _snp_base_codes(snpinfo["snp_base"]);
	IntegerVector snp_base_codes(_snp_base_codes);

	int n_motifs = pwms.size();
	int n_snps = sequence_matrix.ncol();
	
	NumericMatrix log_enhance_odds(n_snps, n_motifs);
	NumericMatrix log_reduce_odds(n_snps, n_motifs);
	IntegerMatrix match_ref_base(n_snps, n_motifs);
	IntegerMatrix match_snp_base(n_snps, n_motifs);
	NumericMatrix log_lik_ratio(n_snps, n_motifs);
	NumericMatrix log_lik_ref(n_snps, n_motifs);
	NumericMatrix log_lik_snp(n_snps, n_motifs);
	NumericMatrix mean_log_lik_ratio(n_snps, n_motifs);
	NumericMatrix mean_log_lik_ref(n_snps, n_motifs);
	NumericMatrix mean_log_lik_snp(n_snps, n_motifs);
	NumericMatrix median_log_lik_ratio(n_snps, n_motifs);
	NumericMatrix median_log_lik_ref(n_snps, n_motifs);
	NumericMatrix median_log_lik_snp(n_snps, n_motifs);
	
	double tol = 1e-10;

	// change all inputs to 0 indexed
	for(int snp_id = 0; snp_id < n_snps; snp_id ++) {
		for(int i = 0; i < sequence_matrix.nrow(); i ++) {
			sequence_matrix(i, snp_id) --;
		}
		ref_base_codes[snp_id] --;
		snp_base_codes[snp_id] --;
	}

	//for each snp
	for(int snp_id = 0; snp_id < n_snps; snp_id ++) {
		//construct reverse sequence
		IntegerVector snp_sequence = sequence_matrix(_, snp_id);
		snp_sequence[snp_sequence.size() / 2] = ref_base_codes[snp_id];
		IntegerVector snp_sequence_snp_base(clone(snp_sequence));
		snp_sequence_snp_base[snp_sequence.size() / 2] = snp_base_codes[snp_id];

		// for each motif
		for(int motif_id = 0; motif_id < n_motifs; motif_id ++) {
			SEXP _pwm(pwms[motif_id]);
			NumericMatrix pwm(_pwm);
			for(int i = 0; i < pwm.nrow(); i ++)
				for(int j = 0; j < pwm.ncol(); j ++)
					if(pwm(i, j) < tol)
						pwm(i, j) = tol;
					
			/*
			int _match_ref_base = find_best_match(pwm, snp_sequence);
			int _match_snp_base = find_best_match(pwm, snp_sequence_snp_base);
			log_lik_ref(snp_id, motif_id) = bidir_pwm_log_prob(pwm, snp_sequence, _match_ref_base);
			log_lik_snp(snp_id, motif_id) = bidir_pwm_log_prob(pwm, snp_sequence_snp_base, _match_snp_base);
			log_lik_ratio(snp_id, motif_id) = log_lik_ref(snp_id, motif_id) - log_lik_snp(snp_id, motif_id);
			match_ref_base(snp_id, motif_id) = _match_ref_base;
			match_snp_base(snp_id, motif_id) = _match_snp_base;
			*/
			SequenceScores ref_scores = comp_seq_scores(pwm, snp_sequence);
			SequenceScores snp_scores = comp_seq_scores(pwm, snp_sequence_snp_base);
			int _match_ref_base = ref_scores.best_match_pos;
			int _match_snp_base = snp_scores.best_match_pos;
			match_ref_base(snp_id, motif_id) = ref_scores.best_match_pos;
			match_snp_base(snp_id, motif_id) = snp_scores.best_match_pos;

      // max log lik
			log_lik_ref(snp_id, motif_id) = ref_scores.max_log_lik;
			log_lik_snp(snp_id, motif_id) = snp_scores.max_log_lik;
			log_lik_ratio(snp_id, motif_id) = ref_scores.max_log_lik - snp_scores.max_log_lik;
			// mean log lik
			mean_log_lik_ref(snp_id, motif_id) = ref_scores.mean_log_lik;
			mean_log_lik_snp(snp_id, motif_id) = snp_scores.mean_log_lik;
			mean_log_lik_ratio(snp_id, motif_id) = ref_scores.mean_log_lik - snp_scores.mean_log_lik;
			// median log lik
			median_log_lik_ref(snp_id, motif_id) = ref_scores.median_log_lik;
			median_log_lik_snp(snp_id, motif_id) = snp_scores.median_log_lik;
			median_log_lik_ratio(snp_id, motif_id) = ref_scores.median_log_lik - snp_scores.median_log_lik;
			
			// computing for ref_base
			if(_match_ref_base > 0) {
			  _match_ref_base --;
				int snp_pos_in_pwm = snp_sequence.size() / 2 - _match_ref_base;
				log_reduce_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, ref_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, snp_base_codes[snp_id]));
				// note: index from 0 is confusing; cannot separate positive from negative
			} else {
			  _match_ref_base = - _match_ref_base - 1;
				int snp_pos_in_pwm = snp_sequence.size() / 2 - _match_ref_base;
				log_reduce_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, 3 - ref_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, 3 - snp_base_codes[snp_id]));
			}
			// computing for snp_base
			if(_match_snp_base > 0) {
			  _match_snp_base --;
				int snp_pos_in_pwm = snp_sequence.size() / 2 - _match_snp_base;
				log_enhance_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, snp_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, ref_base_codes[snp_id]));
			} else {
			  _match_snp_base = - _match_snp_base - 1;
				int snp_pos_in_pwm = snp_sequence.size() / 2 - _match_snp_base;
				log_enhance_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, 3 - snp_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, 3 - ref_base_codes[snp_id])); 
			}
		}
	}

	return Rcpp::List::create(
	  Rcpp::Named("log_enhance_odds") = log_enhance_odds,
	  Rcpp::Named("log_reduce_odds") = log_reduce_odds,
    Rcpp::Named("match_ref_base") = match_ref_base,
    Rcpp::Named("match_snp_base") = match_snp_base,
    Rcpp::Named("log_lik_ratio") = log_lik_ratio,
    Rcpp::Named("log_lik_ref") = log_lik_ref,
    Rcpp::Named("log_lik_snp") = log_lik_snp,
    Rcpp::Named("mean_log_lik_ref") = mean_log_lik_ref,
    Rcpp::Named("mean_log_lik_snp") = mean_log_lik_snp,
    Rcpp::Named("mean_log_lik_ratio") = mean_log_lik_ratio,
    Rcpp::Named("median_log_lik_ref") = median_log_lik_ref,
    Rcpp::Named("median_log_lik_snp") = median_log_lik_snp,
    Rcpp::Named("median_log_lik_ratio") = median_log_lik_ratio
	);
}


int get_min_start_pos(int motif_len, int seq_len) {
  int min_start_pos = seq_len / 2 - motif_len + 1;
  if(min_start_pos < 0)
    min_start_pos = 0;
  return(min_start_pos);
}
  
int get_max_start_pos(int motif_len, int seq_len) {
  int max_start_pos = seq_len / 2;
  if(max_start_pos + motif_len - 1 >= seq_len) {
    max_start_pos = seq_len - motif_len;
  }
  return(max_start_pos);
}

IntegerVector revert_sequence(IntegerVector sequence) {
  int seq_len = sequence.size();
  IntegerVector rev_sequence(seq_len);
  for(int i = 0; i < seq_len; ++ i) {
    rev_sequence[i] = 3 - sequence[seq_len - 1 - i];
  }
  return(rev_sequence);
}

/*
 * Returns a vector of scores for subsequences from both strands.
 * 2*i-th position is the score for the i-th subsequence on the positive strand,
 * (2*i+1)-th position is the score for the i-th subsequence on the reverse strand.
 */
NumericVector comp_subseq_scores(NumericMatrix pwm, IntegerVector sequence) {
  int motif_len = pwm.nrow();
  int seq_len = sequence.size();
  int min_start_pos = get_min_start_pos(motif_len, seq_len);
  int max_start_pos = get_max_start_pos(motif_len, seq_len);
  
  int n_subsequences = max_start_pos - min_start_pos + 1;
  NumericVector subseq_scores(n_subsequences * 2);
  IntegerVector rev_sequence = revert_sequence(sequence);

  for(int start_pos = min_start_pos; start_pos <= max_start_pos; start_pos ++) {
    int idx = 2 * (start_pos - min_start_pos);
    subseq_scores[idx] = pwm_log_prob(pwm, sequence, start_pos);
    subseq_scores[idx + 1] = pwm_log_prob(pwm, rev_sequence, start_pos);
  }
  return subseq_scores;
}

/*
 * Find the subsequence in a "sequence" that best matches the position weight matrix ("pwm").
 * If the best subsequence has starting position x on the positive strand, this returns x+1;
 * If the best subsequence has starting position x on the negative strand, this returns -x-1.
 * Returns the struct containing (best match position, max loglik, mean loglik, median loglik).
 */
SequenceScores comp_seq_scores(NumericMatrix pwm, IntegerVector sequence) {
  NumericVector subseq_scores = comp_subseq_scores(pwm, sequence);
  int raw_max_idx = Rcpp::which_max(subseq_scores);
  int match_pos = raw_max_idx / 2 + get_min_start_pos(pwm.nrow(), sequence.size());
  match_pos ++;
  if(raw_max_idx % 2 == 1) {
    // This is the negative strand
    match_pos = -match_pos;
  }
  float max_log_lik = Rcpp::max(subseq_scores);
  float mean_log_lik = Rcpp::mean(subseq_scores);
  float median_log_lik = Rcpp::mean(subseq_scores);
  SequenceScores scores = {
    match_pos, max_log_lik, mean_log_lik, median_log_lik
  };
  return scores;
}


/* 
 * Compute the log likelihood corresponding to a position weight matrix ("pwm") for 
 * a subsequence in "sequence" starting from "start_pos".
*/
double pwm_log_prob(NumericMatrix pwm, IntegerVector sequence, int start_pos) {
  assert(start_pos >= 0);
	double tol = 1e-10;
	double log_prob = 0;
	for(int i = 0; i < pwm.nrow(); i ++)
		for(int j = 0; j < pwm.ncol(); j ++)
			if(pwm(i, j) < tol)
				pwm(i, j) = tol;
	for(int pos = start_pos; pos < start_pos + pwm.nrow(); pos ++) {
		log_prob += log(pwm(pos - start_pos, sequence[pos]));
	}
	return(log_prob);
}

/* 
 * Compute the log likelihood corresponding to a position weight matrix ("pwm") for 
 * a subsequence in "sequence".
 * If start_pos > 0, this computes for the subsequence starting from start_pos-1 on the positive strand.
 * If start_pos < 0, this computes for the subsequence starting from -start_pos-1 on the reverse strand.
 */
double bidir_pwm_log_prob(NumericMatrix pwm, IntegerVector sequence, int start_pos) {
  assert(start_pos != 0);
  IntegerVector _sequence(sequence);
  if(start_pos < 0) {
    _sequence = revert_sequence(sequence);
    start_pos = -start_pos;
  }
  return(pwm_log_prob(pwm, _sequence, start_pos-1));
}


SEXP test_max_log_prob(SEXP _pwm, SEXP _sequence) {
	NumericMatrix pwm(_pwm);
	IntegerVector sequence(_sequence);
	return(wrap(comp_seq_scores(pwm, sequence).max_log_lik));
}

/*
Compute the first order Markovian transition matrix.
*/
SEXP transition_matrix(SEXP _sequence_matrix) {
	IntegerMatrix sequence_matrix(_sequence_matrix);
	NumericMatrix transition_matrix(4, 4);
	for(int i = 0; i < sequence_matrix.nrow() - 1; i ++) {
		for(int j = 0; j < sequence_matrix.ncol(); j ++) {
			transition_matrix(sequence_matrix(i, j) - 1, sequence_matrix(i + 1, j) - 1) ++;
		}
	}
	return(wrap(transition_matrix));
}

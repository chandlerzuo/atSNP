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
SEXP motif_score(SEXP _motif_library, SEXP _snpinfo) {
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
	NumericMatrix log_lik_ratio(n_snps, n_motifs);
	IntegerMatrix match_ref_base(n_snps, n_motifs);
	IntegerMatrix match_snp_base(n_snps, n_motifs);
	NumericMatrix log_lik_ref(n_snps, n_motifs);
	NumericMatrix log_lik_snp(n_snps, n_motifs);

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
		//printf("snpid: %d\n", snp_id);
		//construct reverse sequence
		IntegerVector snp_sequence = sequence_matrix(_, snp_id);
		snp_sequence[snp_sequence.size() / 2] = ref_base_codes[snp_id];
		IntegerVector rev_sequence(sequence_matrix.nrow());
		for(int i = 0; i < rev_sequence.size(); i ++) {
			rev_sequence[i] = 3 - snp_sequence[snp_sequence.size() - 1 - i];
		}
		IntegerVector snp_sequence_snp_base(clone(snp_sequence));
		IntegerVector rev_sequence_snp_base(clone(rev_sequence));
		snp_sequence_snp_base[snp_sequence.size() / 2] = snp_base_codes[snp_id];
		rev_sequence_snp_base[snp_sequence.size() / 2] = 3 - snp_base_codes[snp_id];
		// for each motif
		for(int motif_id = 0; motif_id < n_motifs; motif_id ++) {
			SEXP _pwm(pwms[motif_id]);
			NumericMatrix pwm(_pwm);
			for(int i = 0; i < pwm.nrow(); i ++)
				for(int j = 0; j < pwm.ncol(); j ++)
					if(pwm(i, j) < tol)
						pwm(i, j) = tol;
			// find maximum on the positive strand, ref_base
			int match_pos_ref_base = find_best_match(pwm, snp_sequence);
			double log_prob_pos_ref_base = pwm_log_prob(pwm, snp_sequence, match_pos_ref_base);
			// find maximum on the reverse strand, ref_base
			int match_rev_ref_base = find_best_match(pwm, rev_sequence);
			double log_prob_rev_ref_base = pwm_log_prob(pwm, rev_sequence, match_rev_ref_base);
			// find maximum on the positive strand, snp_base
			int match_pos_snp_base = find_best_match(pwm, snp_sequence_snp_base);
			double log_prob_pos_snp_base = pwm_log_prob(pwm, snp_sequence_snp_base, match_pos_snp_base);
			// find maximum on the reverse strand, snp_base
			int match_rev_snp_base = find_best_match(pwm, rev_sequence_snp_base);
			double log_prob_rev_snp_base = pwm_log_prob(pwm, rev_sequence_snp_base, match_rev_snp_base);
			// computing for ref_base
			double log_prob_ref = 0;
			if(log_prob_pos_ref_base > log_prob_rev_ref_base) {
				log_prob_ref = log_prob_pos_ref_base;
				int snp_pos_in_pwm = snp_sequence.size() / 2 - match_pos_ref_base;
				log_reduce_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, ref_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, snp_base_codes[snp_id]));
				// note: index from 0 is confusing; cannot separate positive from negative
				match_ref_base(snp_id, motif_id) = match_pos_ref_base + 1;
			} else {
				log_prob_ref = log_prob_rev_ref_base;
				int snp_pos_in_pwm = rev_sequence.size() / 2 - match_rev_ref_base;
				log_reduce_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, 3 - ref_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, 3 - snp_base_codes[snp_id])); 
				match_ref_base(snp_id, motif_id) = - match_rev_ref_base - 1;
			}
			log_lik_ref(snp_id, motif_id) = log_prob_ref;
			// computing for snp_base
			double log_prob_snp = 0;
			if(log_prob_pos_snp_base > log_prob_rev_snp_base) {
				log_prob_snp = log_prob_pos_snp_base;
				int snp_pos_in_pwm = snp_sequence.size() / 2 - match_pos_snp_base;
				log_enhance_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, snp_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, ref_base_codes[snp_id])); 
				match_snp_base(snp_id, motif_id) = match_pos_snp_base + 1;
			} else {
				log_prob_snp = log_prob_rev_snp_base;
				int snp_pos_in_pwm = rev_sequence.size() / 2 - match_rev_snp_base;
				log_enhance_odds(snp_id, motif_id) = log(pwm(snp_pos_in_pwm, 3 - snp_base_codes[snp_id])) - 
					log(pwm(snp_pos_in_pwm, 3 - ref_base_codes[snp_id])); 
				match_snp_base(snp_id, motif_id) = - match_rev_snp_base - 1;
			}
			log_lik_snp(snp_id, motif_id) = log_prob_snp;
			// log likelihood ratio
			log_lik_ratio(snp_id, motif_id) = log_prob_ref - log_prob_snp;
		}
	}

	return Rcpp::List::create(
				  Rcpp::Named("log_enhance_odds") = log_enhance_odds,
				  Rcpp::Named("log_reduce_odds") = log_reduce_odds,
				  Rcpp::Named("match_ref_base") = match_ref_base,
				  Rcpp::Named("match_snp_base") = match_snp_base,
				  Rcpp::Named("log_lik_ratio") = log_lik_ratio,
				  Rcpp::Named("log_lik_ref") = log_lik_ref,
				  Rcpp::Named("log_lik_snp") = log_lik_snp);
}

/*
Find the subsequence in a "sequence" that best matches the position weight matrix ("pwm"). Returns the starting position of this subsequence.
*/
int find_best_match(NumericMatrix pwm, IntegerVector sequence) {
	int motif_len = pwm.nrow();
	int seq_len = sequence.size();

	int min_start_pos = seq_len / 2 - motif_len + 1;
	if(min_start_pos < 0)
		min_start_pos = 0;

	int max_start_pos = seq_len / 2;
	if(max_start_pos + motif_len - 1 >= seq_len) {
		max_start_pos = seq_len - motif_len;
	}

	double max_log_prob = -100 * motif_len;
	double log_prob = max_log_prob;
	int max_pos = min_start_pos;
	for(int start_pos = min_start_pos; start_pos <= max_start_pos; start_pos ++) {
		log_prob = pwm_log_prob(pwm, sequence, start_pos);
		if(log_prob > max_log_prob) {
			max_log_prob = log_prob;
			max_pos = start_pos;
		}
		// printf("start %d, log_prob %lf, max so far %lf\n", start_pos, log_prob, max_log_prob);
	}
	return(max_pos);
}

/* 
Compute the log likelihood corresponding to a position weight matrix ("pwm") for a subsequence in "sequence" starting from "start_pos".
*/
double pwm_log_prob(NumericMatrix pwm, IntegerVector sequence, int start_pos) {
	double tol = 1e-10;
	double log_prob = 0;
	for(int i = 0; i < pwm.nrow(); i ++)
		for(int j = 0; j < pwm.ncol(); j ++)
			if(pwm(i, j) < tol)
				pwm(i, j) = tol;
	for(int pos = start_pos; pos < start_pos + pwm.nrow(); pos ++) {
		//printf("pwm(%d-%d,%d)=%3.3f\n", pos, start_pos, sequence[pos] - 1, pwm(pos - start_pos, sequence[pos] - 1));
		log_prob += log(pwm(pos - start_pos, sequence[pos]) );
	}
	return(log_prob);
}

SEXP test_max_log_prob(SEXP _pwm, SEXP _sequence) {
	NumericMatrix pwm(_pwm);
	IntegerVector sequence(_sequence);
	return(wrap(pwm_log_prob(pwm, sequence, find_best_match(pwm, sequence))));
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

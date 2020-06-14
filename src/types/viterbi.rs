use super::constants::*;
use super::train::{get_protein, get_rc_dna, get_rc_dna_indel, nt2int, trinucleotide, Train, HMM};
use std::convert::TryInto;
use std::fmt::Write;

#[derive(Clone, Debug)]
pub struct Output {
	pub out: String,
	pub aa: String,
	pub dna: String
}

pub fn viterbi(
	hmm: &HMM,
	train: &Train,
	O: &String,
	wholegenome: bool,
	cg: usize,
	format: bool,
	head: &String,
) -> Output {
	let log53: f64 = (0.53 as f64).ln();
	let log16: f64 = (0.16 as f64).ln();
	let log30: f64 = (0.30 as f64).ln();
	let log25: f64 = (0.25 as f64).ln();
	let log95: f64 = (0.95 as f64).ln();
	let log54: f64 = (0.54 as f64).ln();
	let log83: f64 = (0.83 as f64).ln();
	let log07: f64 = (0.07 as f64).ln();
	let max_dbl = f64::INFINITY;

	let gene_len;

	let mut dna_id = 0;
	let mut dna_f_id = 0;

	let mut temp_i = [0; 6];
	let mut temp_i_1 = [0; 6];

	let mut refine = false;

	if wholegenome {
		gene_len = 120;
		refine = true;
	} else {
		gene_len = 60;
	}

	let len_seq = O.chars().count();
	let mut alpha = vec![vec![0.0; len_seq]; NUM_STATE];
	let mut path: Vec<Vec<i8>> = vec![vec![0; len_seq]; NUM_STATE];
	let mut vpath = vec![0; len_seq];

	for i in 0..NUM_STATE {
		alpha[i][0] = -1.0 * hmm.initial_state[i];
	}

	let O: &Vec<char> = &O.chars().collect();

	/* stop state */
	if (O[0] == 'T' || O[0] == 't')
		&& (((O[1] == 'A' || O[1] == 'a') && (O[2] == 'A' || O[2] == 'a'))
			|| ((O[1] == 'A' || O[1] == 'a') && (O[2] == 'G' || O[2] == 'g'))
			|| ((O[1] == 'G' || O[1] == 'g') && (O[2] == 'A' || O[2] == 'a')))
	// TODO: this needs to be a lot cleaner!
	{
		alpha[E_STATE][0] = max_dbl;
		alpha[E_STATE][1] = max_dbl;
		path[E_STATE][1] = E_STATE as i8;
		path[E_STATE][2] = E_STATE as i8;

		alpha[M6_STATE][2] = max_dbl;
		alpha[M5_STATE][1] = max_dbl;
		alpha[M4_STATE][0] = max_dbl;
		alpha[M3_STATE][2] = max_dbl;
		alpha[M2_STATE][1] = max_dbl;
		alpha[M1_STATE][0] = max_dbl;

		if (O[1] == 'A' || O[1] == 'a') && (O[2] == 'A' || O[2] == 'a') {
			alpha[E_STATE][2] = alpha[E_STATE][2] - log53;
		} else if (O[1] == 'A' || O[1] == 'a') && (O[2] == 'G' || O[2] == 'g') {
			alpha[E_STATE][2] = alpha[E_STATE][2] - log16;
		} else if (O[1] == 'G' || O[1] == 'g') && (O[2] == 'A' || O[2] == 'a') {
			alpha[E_STATE][2] = alpha[E_STATE][2] - log30;
		}
	}

	if (O[2] == 'A' || O[2] == 'a')
		&& (((O[0] == 'T' || O[0] == 't') && (O[1] == 'T' || O[1] == 't'))
			|| ((O[0] == 'C' || O[0] == 'c') && (O[1] == 'T' || O[1] == 't'))
			|| ((O[0] == 'T' || O[0] == 't') && (O[1] == 'C' || O[1] == 'c')))
	{
		alpha[S_STATE_1][0] = max_dbl;
		alpha[S_STATE_1][1] = max_dbl;
		alpha[S_STATE_1][2] = alpha[S_STATE][0];
		path[S_STATE_1][1] = S_STATE_1 as i8;
		path[S_STATE_1][2] = S_STATE_1 as i8;

		alpha[M3_STATE_1][2] = max_dbl;
		alpha[M6_STATE_1][2] = max_dbl;

		if (O[0] == 'T' || O[0] == 't') && (O[1] == 'T' || O[1] == 't') {
			alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log53;
		} else if (O[0] == 'C' || O[0] == 'c') && (O[1] == 'T' || O[1] == 't') {
			alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log16;
		} else if (O[0] == 'T' || O[0] == 't') && (O[1] == 'C' || O[1] == 'c') {
			alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log30;
		}
	}
	/******************************************************************/
	/*  fill out the rest of the columns                              */
	/******************************************************************/
	let mut num_N = 0;
	for t in 1..len_seq {
		let mut from = nt2int(O[t - 1]);
		let mut from0: usize;
		if t > 2 {
			from0 = nt2int(O[t - 2]);
		}
		else{
			from0 = 2;
		}
		let mut to = nt2int(O[t]);

		/* if DNA is other than ACGT, do it later */
		if from == 4 {
			from = 2;
		}
		if from0 == 4 {
			from0 = 2;
		}
		if to == 4 {
			to = 2;
			num_N += 1;
		} else {
			num_N = 0;
		}
		let from2 = from0 * 4 + from;

		/******************/
		/* M state        */
		/******************/

		for i in M1_STATE..=M6_STATE {
			if alpha[i][t] < max_dbl {
				if t != 0 {
					let mut j;
					if i == M1_STATE {
						/* from M state */
						j = M6_STATE;
						alpha[i][t] =
							alpha[j][t - 1] - hmm.tr[TR_GG] - hmm.tr[TR_MM] - hmm.e_m[0][from2][to];
						path[i][t] = j as i8;

						/* from D state */
						if !wholegenome {
							for j in (M1_STATE..=M5_STATE).rev() {
								let mut num_d: isize = -10;
								if j >= i {
									num_d = i as isize - j as isize + 6; // TODO: this seems ugly
								} else if j + 1 < i {
									num_d = i as isize - j as isize; // TODO: this seems ugly
								}
								if num_d > 0 {
									let temp_alpha = alpha[j][t - 1]
										- hmm.tr[TR_MD] - hmm.e_m[0][from2][to] - log25
										* (num_d as f64 - 1.0) - hmm.tr[TR_DD]
										* (num_d as f64 - 2.0) - hmm.tr[TR_DM];
									if temp_alpha < alpha[i][t] {
										alpha[i][t] = temp_alpha;
										path[i][t] = j as i8;
									}
								}
							}
						}

						/* from Start state */
						let temp_alpha = alpha[S_STATE][t - 1] - hmm.e_m[0][from2][to];
						if temp_alpha < alpha[i][t] {
							alpha[i][t] = temp_alpha;
							path[i][t] = S_STATE as i8;
						}
					} else {
						/*i ==M2-M6*/

						/* from M state */
						j = i - 1;
						alpha[i][t] =
							alpha[j][t - 1] - hmm.tr[TR_MM] - hmm.e_m[i - M1_STATE][from2][to];
						path[i][t] = j as i8;

						/* from D state */
						if !wholegenome {
							for j in (M1_STATE..=M6_STATE).rev() {
								let mut num_d: isize = -10;
								if j >= i {
									num_d = i as isize - j as isize + 6; // TODO: this seems ugly
								} else if j + 1 < i {
									num_d = i as isize - j as isize; // TODO: this seems ugly
								}
								if num_d > 0 {
									let temp_alpha = alpha[j][t - 1]
										- hmm.tr[TR_MD] - hmm.e_m[i - M1_STATE][from2]
										[to] - log25 * (num_d as f64 - 1.0) - hmm.tr
										[TR_DD]
										* (num_d as f64 - 2.0) - hmm.tr[TR_DM];
									if temp_alpha < alpha[i][t] {
										alpha[i][t] = temp_alpha;
										path[i][t] = j as i8;
									}
								}
							}
						}
					}

					/* from I state */
					if i == M1_STATE {
						j = I6_STATE;
					} else {
						j = I1_STATE + (i - M1_STATE - 1);
					}

					/* to avoid stop codon */
					if t < 2 {
					} else if (i == M2_STATE || i == M5_STATE)
						&& (O[temp_i[j - I1_STATE]] == 'T' || O[temp_i[j - I1_STATE]] == 't')
						&& t < len_seq - 1 && ( ((O[t] == 'A' || O[t] == 'a') && (O[t + 1] == 'A' || O[t + 1] == 'a'))
							|| ((O[t] == 'A' || O[t] == 'a')
								&& (O[t + 1] == 'G' || O[t + 1] == 'g'))
							|| ((O[t] == 'G' || O[t] == 'g')
								&& (O[t + 1] == 'A' || O[t + 1] == 'a')))
					{
					}
					else if (i == M3_STATE || i == M6_STATE) && temp_i[j - I1_STATE] as isize - 1 > 0
						&& (O[temp_i[j - I1_STATE] - 1] == 'T'
							|| O[temp_i[j - I1_STATE] - 1] == 't')
						&& (((O[temp_i[j - I1_STATE]] == 'A' || O[temp_i[j - I1_STATE]] == 'a')
							&& (O[t] == 'A' || O[t] == 'a'))
							|| ((O[temp_i[j - I1_STATE]] == 'A' || O[temp_i[j - I1_STATE]] == 'a')
								&& (O[t] == 'G' || O[t] == 'g'))
							|| ((O[temp_i[j - I1_STATE]] == 'G' || O[temp_i[j - I1_STATE]] == 'g')
								&& (O[t] == 'A' || O[t] == 'a')))
					{
					} else {
						let temp_alpha = alpha[j][t - 1] - hmm.tr[TR_IM] - log25;
						if temp_alpha < alpha[i][t] {
							alpha[i][t] = temp_alpha;
							path[i][t] = j as i8;
						}
					}
				}
			}
		}

		/******************/
		/* I state        */
		/******************/
		for i in I1_STATE..=I6_STATE {
			let mut j;
			if t != 0 {
				/* from I state */
				j = i;
				alpha[i][t] = alpha[j][t - 1] - hmm.tr[TR_II] - hmm.tr_i_i[from][to];
				path[i][t] = j as i8;

				/* from M state */
				j = i - I1_STATE + M1_STATE;
				let temp_alpha;
				if i == I6_STATE {
					temp_alpha =
						alpha[j][t - 1] - hmm.tr[TR_GG] - hmm.tr[TR_MI] - hmm.tr_m_i[from][to];
				} else {
					temp_alpha = alpha[j][t - 1] - hmm.tr[TR_MI] - hmm.tr_m_i[from][to];
				}
				if temp_alpha < alpha[i][t] {
					alpha[i][t] = temp_alpha;
					path[i][t] = j as i8;

					temp_i[i - I1_STATE] = t - 1;
				}
			}
		}

		/******************/
		/* M' state        */
		/******************/

		for i in M1_STATE_1..=M6_STATE_1 {
			let mut j;
			if (i == M1_STATE_1 || i == M4_STATE_1)
				&& t >= 3 && (((O[t - 3] == 'T' || O[t - 3] == 't')
				&& (O[t - 2] == 'T' || O[t - 2] == 't')
				&& (O[t - 1] == 'A' || O[t - 1] == 'a'))
				|| ((O[t - 3] == 'C' || O[t - 3] == 'c')
					&& (O[t - 2] == 'T' || O[t - 2] == 't')
					&& (O[t - 1] == 'A' || O[t - 1] == 'a'))
				|| ((O[t - 3] == 'T' || O[t - 3] == 't')
					&& (O[t - 2] == 'C' || O[t - 2] == 'c')
					&& (O[t - 1] == 'A' || O[t - 1] == 'a')))
			{
				/* from Start state  since this is actually stop codon in minus strand */
				alpha[i][t] = alpha[S_STATE_1][t - 1] - hmm.e_m_1[i - M1_STATE_1][from2][to];
				path[i][t] = S_STATE_1 as i8;
			} else {
				if t != 0 {
					if i == M1_STATE_1 {
						/* from M state */
						j = M6_STATE_1;
						alpha[i][t] = alpha[j][t - 1]
							- hmm.tr[TR_GG] - hmm.tr[TR_MM]
							- hmm.e_m_1[0][from2][to];
						path[i][t] = j as i8;

						/* from D state */
						if !wholegenome {
							for j in (M1_STATE_1..=M5_STATE_1).rev() {
								let mut num_d: isize = -10;
								if j >= i {
									num_d = i as isize - j as isize + 6; // TODO: this seems ugly
								} else if j + 1 < i {
									num_d = i as isize - j as isize; // TODO: this seems ugly
								}
								if num_d > 0 {
									let temp_alpha = alpha[j][t - 1]
										- hmm.tr[TR_MD] - hmm.e_m_1[0][from2][to]
										- log25 * (num_d as f64 - 1.0) - hmm.tr[TR_DD]
										* (num_d as f64 - 2.0) - hmm.tr[TR_DM];
									if temp_alpha < alpha[i][t] {
										alpha[i][t] = temp_alpha;
										path[i][t] = j as i8;
									}
								}
							}
						}
					} else {
						/* from M state */
						j = i - 1;
						alpha[i][t] =
							alpha[j][t - 1] - hmm.tr[TR_MM] - hmm.e_m_1[i - M1_STATE_1][from2][to];
						path[i][t] = j as i8;

						/* from D state */
						if !wholegenome {
							for j in (M1_STATE_1..=M6_STATE_1).rev() {
								let mut num_d: isize = -10;
								if j >= i {
									num_d = i as isize - j as isize + 6; // TODO: this seems ugly
								} else if j + 1 < i {
									num_d = i as isize - j as isize; // TODO: this seems ugly
								}
								if num_d > 0 {
									let temp_alpha = alpha[j][t - 1]
										- hmm.tr[TR_MD] - hmm.e_m_1[i - M1_STATE_1]
										[from2][to] - log25 * (num_d as f64 - 1.0)
										- hmm.tr[TR_DD] * (num_d as f64 - 2.0) - hmm
										.tr[TR_DM];
									if temp_alpha < alpha[i][t] {
										alpha[i][t] = temp_alpha;
										path[i][t] = j as i8;
									}
								}
							}
						}
					}

					/* from I state */
					if i == M1_STATE_1 {
						j = I6_STATE_1;
					} else {
						j = I1_STATE_1 + (i - M1_STATE_1 - 1);
					}

					/* to avoid stop codon */
					if t < 2 || t == len_seq - 1 {
					} else if (i == M2_STATE_1 || i == M5_STATE_1)
						&& (O[t + 1] == 'A' || O[t + 1] == 'a')
						&& (((O[temp_i_1[j - I1_STATE_1]] == 'T'
							|| O[temp_i_1[j - I1_STATE_1]] == 't')
							&& (O[t] == 'T' || O[t] == 't'))
							|| ((O[temp_i_1[j - I1_STATE_1]] == 'C'
								|| O[temp_i_1[j - I1_STATE_1]] == 'c')
								&& (O[t] == 'T' || O[t] == 't'))
							|| ((O[temp_i_1[j - I1_STATE_1]] == 'T'
								|| O[temp_i_1[j - I1_STATE_1]] == 't')
								&& (O[t] == 'C' || O[t] == 'c')))
					{
					} else if (i == M3_STATE_1 || i == M6_STATE_1)
						&& (O[t] == 'A' || O[t] == 'a') && temp_i_1[j - I1_STATE_1] as isize - 1 > 0
						&& (((O[temp_i_1[j - I1_STATE_1] - 1] == 'T'
							|| O[temp_i_1[j - I1_STATE_1] - 1] == 't')
							&& (O[temp_i_1[j - I1_STATE_1]] == 'T'
								|| O[temp_i_1[j - I1_STATE_1]] == 't'))
							|| ((O[temp_i_1[j - I1_STATE_1] - 1] == 'C'
								|| O[temp_i_1[j - I1_STATE_1] - 1] == 'c')
								&& (O[temp_i_1[j - I1_STATE_1]] == 'T'
									|| O[temp_i_1[j - I1_STATE_1]] == 't'))
							|| ((O[temp_i_1[j - I1_STATE_1] - 1] == 'T'
								|| O[temp_i_1[j - I1_STATE_1] - 1] == 't')
								&& (O[temp_i_1[j - I1_STATE_1]] == 'C'
									|| O[temp_i_1[j - I1_STATE_1]] == 'c')))
					{
					} else {
						let temp_alpha = alpha[j][t - 1] - hmm.tr[TR_IM] - log25;
						if temp_alpha < alpha[i][t] {
							alpha[i][t] = temp_alpha;
							path[i][t] = j as i8;
						}
					}
				}
			}
		}

		/******************/
		/* I' state        */
		/******************/
		for i in I1_STATE_1..=I6_STATE_1 {
			if t != 0 {
				/* from I state */
				let mut j = i;
				alpha[i][t] = alpha[j][t - 1] - hmm.tr[TR_II] - hmm.tr_i_i[from][to];
				path[i][t] = j as i8;

				/* from M state */
				if t > 4 && path[S_STATE_1][t - 3] != R_STATE as i8
					&& path[S_STATE_1][t - 4] != R_STATE as i8
					&& path[S_STATE_1][t - 5] != R_STATE as i8
				{
					j = i - I1_STATE_1 + M1_STATE_1;
					let temp_alpha;
					if i == I6_STATE_1 {
						temp_alpha =
							alpha[j][t - 1] - hmm.tr[TR_GG] - hmm.tr[TR_MI] - hmm.tr_m_i[from][to];
					} else {
						temp_alpha = alpha[j][t - 1] - hmm.tr[TR_MI] - hmm.tr_m_i[from][to];
					}
					if temp_alpha < alpha[i][t] {
						alpha[i][t] = temp_alpha;
						path[i][t] = j as i8;

						temp_i_1[i - I1_STATE_1] = t - 1;
					}
				}
			}
		}

		/***********************/
		/* Non_coding state    */
		/***********************/

		if t != 0 {
			alpha[R_STATE][t] = alpha[R_STATE][t - 1] - hmm.tr_r_r[from][to] - hmm.tr[TR_RR];
			path[R_STATE][t] = R_STATE as i8;

			let mut temp_alpha = alpha[E_STATE][t - 1] - hmm.tr[TR_ER];
			if temp_alpha < alpha[R_STATE][t] {
				alpha[R_STATE][t] = temp_alpha;
				path[R_STATE][t] = E_STATE as i8;
			}

			temp_alpha = alpha[E_STATE_1][t - 1] - hmm.tr[TR_ER];
			if temp_alpha < alpha[R_STATE][t] {
				alpha[R_STATE][t] = temp_alpha;
				path[R_STATE][t] = E_STATE_1 as i8;
			}
			alpha[R_STATE][t] -= log95;
		}

		/******************/
		/* END state      */
		/******************/
		if alpha[E_STATE][t] == 0.0 {
			alpha[E_STATE][t] = max_dbl;
			path[E_STATE][t] = NOSTATE;

			if t < len_seq - 2
				&& (O[t] == 'T' || O[t] == 't')
				&& (((O[t + 1] == 'A' || O[t + 1] == 'a') && (O[t + 2] == 'A' || O[t + 2] == 'a'))
					|| ((O[t + 1] == 'A' || O[t + 1] == 'a')
						&& (O[t + 2] == 'G' || O[t + 2] == 'g'))
					|| ((O[t + 1] == 'G' || O[t + 1] == 'g')
						&& (O[t + 2] == 'A' || O[t + 2] == 'a')))
			{
				alpha[E_STATE][t + 2] = max_dbl;
				/* transition from frame4,frame5,and frame6 */
				let mut temp_alpha = alpha[M6_STATE][t - 1] - hmm.tr[TR_GE];
				if temp_alpha < alpha[E_STATE][t + 2] {
					alpha[E_STATE][t + 2] = temp_alpha;
					path[E_STATE][t] = M6_STATE as i8;
				}

				/* transition from frame1,frame2,and frame3 */
				temp_alpha = alpha[M3_STATE][t - 1] - hmm.tr[TR_GE];
				if temp_alpha < alpha[E_STATE][t + 2] {
					alpha[E_STATE][t + 2] = temp_alpha;
					path[E_STATE][t] = M3_STATE as i8;
				}

				alpha[E_STATE][t] = max_dbl;
				alpha[E_STATE][t + 1] = max_dbl;
				path[E_STATE][t + 1] = E_STATE as i8;
				path[E_STATE][t + 2] = E_STATE as i8;

				alpha[M6_STATE][t + 2] = max_dbl;
				alpha[M5_STATE][t + 1] = max_dbl;
				alpha[M4_STATE][t] = max_dbl;
				alpha[M3_STATE][t + 2] = max_dbl;
				alpha[M2_STATE][t + 1] = max_dbl;
				alpha[M1_STATE][t] = max_dbl;

				if (O[t + 1] == 'A' || O[t + 1] == 'a') && (O[t + 2] == 'A' || O[t + 2] == 'a') {
					alpha[E_STATE][t + 2] = alpha[E_STATE][t + 2] - log54;
				} else if (O[t + 1] == 'A' || O[t + 1] == 'a')
					&& (O[t + 2] == 'G' || O[t + 2] == 'g')
				{
					alpha[E_STATE][t + 2] = alpha[E_STATE][t + 2] - log16;
				} else if (O[t + 1] == 'G' || O[t + 1] == 'g')
					&& (O[t + 2] == 'A' || O[t + 2] == 'a')
				{
					alpha[E_STATE][t + 2] = alpha[E_STATE][t + 2] - log30;
				}

				/* adjustment based on probability distribution */
				let mut start_freq = 0.0;

				let mut sub_sum = 0.0;

				if t >= 60 {
					/* bug reported by Yu-Wei */
					for i in (3..=60).rev() {
						let i = -(i as isize);
						if (t as isize) + i + 2 < len_seq as isize
						// TODO: probably needs some more work
						{
							let idx: usize = (i + 60) as usize;
							start_freq -= hmm.tr_e[idx][trinucleotide(
								&O[(t as isize + i) as usize],
								&O[(t as isize + i + 1) as usize],
								&O[(t as isize + i + 2) as usize],
							)];
						}
					}
				} else {
					for i in (3..=t).rev() {
						let i = -(i as isize);
						if t as isize + i + 2 < len_seq as isize
						// TODO: probably needs some more work
						{
							let idx: usize = (i + 60) as usize;
							sub_sum += hmm.tr_e[idx][trinucleotide(
								&O[(t as isize + i) as usize],
								&O[(t as isize + i + 1) as usize],
								&O[(t as isize + i + 2) as usize],
							)];
						}
					}
					sub_sum = sub_sum * 58.0 / (-3.0 + t as f64 + 1.0);
					start_freq -= sub_sum;
				}

				let h_kd = hmm.e_dist[2]
					* (-1.0 * (start_freq - hmm.e_dist[1]).powi(2) / (2.0 * hmm.e_dist[0].powi(2)))
						.exp();
				let r_kd = hmm.e_dist[5]
					* (-1.0 * (start_freq - hmm.e_dist[4]).powi(2) / (2.0 * hmm.e_dist[3].powi(2)))
						.exp();
				let mut p_kd = h_kd / (h_kd + r_kd);
				if p_kd < 0.01 {
					p_kd = 0.01;
				} else if p_kd > 0.99 {
					p_kd = 0.99;
				}
				alpha[E_STATE][t + 2] = alpha[E_STATE][t + 2] - p_kd.ln();
			}
		}

		/*************************************************/
		/* START' state                                  */
		/* origianlly stop codon of genes in - strand    */
		/*************************************************/
		if alpha[S_STATE_1][t] == 0.0 {
			alpha[S_STATE_1][t] = max_dbl;
			path[S_STATE_1][t] = NOSTATE;

			if t < len_seq - 2
				&& (O[t + 2] == 'A' || O[t + 2] == 'a')
				&& (((O[t] == 'T' || O[t] == 't') && (O[t + 1] == 'T' || O[t + 1] == 't'))
					|| ((O[t] == 'C' || O[t] == 'c') && (O[t + 1] == 'T' || O[t + 1] == 't'))
					|| ((O[t] == 'T' || O[t] == 't') && (O[t + 1] == 'C' || O[t + 1] == 'c')))
			{
				alpha[S_STATE_1][t] = max_dbl;
				path[S_STATE_1][t] = R_STATE as i8;
				alpha[S_STATE_1][t + 1] = max_dbl;
				alpha[S_STATE_1][t + 2] = alpha[R_STATE][t - 1] - hmm.tr[TR_RS];
				path[S_STATE_1][t + 1] = S_STATE_1 as i8;
				path[S_STATE_1][t + 2] = S_STATE_1 as i8;

				let mut temp_alpha = alpha[E_STATE_1][t - 1] - hmm.tr[TR_ES];
				if temp_alpha < alpha[S_STATE_1][t + 2] {
					alpha[S_STATE_1][t + 2] = temp_alpha;
					path[S_STATE_1][t] = E_STATE_1 as i8;
				}

				temp_alpha = alpha[E_STATE][t - 1] - hmm.tr[TR_ES1];
				if temp_alpha < alpha[S_STATE_1][t + 2] {
					alpha[S_STATE_1][t + 2] = temp_alpha;
					path[S_STATE_1][t] = E_STATE as i8;
				}

				alpha[M3_STATE_1][t + 2] = max_dbl;
				alpha[M6_STATE_1][t + 2] = max_dbl;

				if (O[t] == 'T' || O[t] == 't') && (O[t + 1] == 'T' || O[t + 1] == 't') {
					alpha[S_STATE_1][t + 2] = alpha[S_STATE_1][t + 2] - log54;
				} else if (O[t] == 'C' || O[t] == 'c') && (O[t + 1] == 'T' || O[t + 1] == 't') {
					alpha[S_STATE_1][t + 2] = alpha[S_STATE_1][t + 2] - log16;
				} else if (O[t] == 'T' || O[t] == 't') && (O[t + 1] == 'C' || O[t + 1] == 'c') {
					alpha[S_STATE_1][t + 2] = alpha[S_STATE_1][t + 2] - log30;
				}

				/* adjustment based on probability distribution */
				let mut start_freq = 0.0;
				for i in 3..=60 {
					if t + i + 2 < len_seq {
						start_freq -= hmm.tr_s_1[i - 3]
							[trinucleotide(&O[t + i], &O[t + i + 1], &O[t + i + 2])];
					}
				}
				let h_kd = hmm.s1_dist[2]
					* (-1.0 * (start_freq - hmm.s1_dist[1]).powi(2)
						/ (2.0 * hmm.s1_dist[0].powi(2)))
					.exp();
				let r_kd = hmm.s1_dist[5]
					* (-1.0 * (start_freq - hmm.s1_dist[4]).powi(2)
						/ (2.0 * hmm.s1_dist[3].powi(2)))
					.exp();
				let mut p_kd = h_kd / (h_kd + r_kd);
				if p_kd < 0.01 {
					p_kd = 0.01;
				} else if p_kd > 0.99 {
					p_kd = 0.99;
				}
				alpha[S_STATE_1][t + 2] = alpha[S_STATE_1][t + 2] - (p_kd).ln();
			}
		}

		/************************/
		/* START state          */
		/************************/
		if alpha[S_STATE][t] == 0.0 {
			alpha[S_STATE][t] = max_dbl;
			path[S_STATE][t] = NOSTATE;

			if t < len_seq - 2
				&& (O[t + 1] == 'T' || O[t + 1] == 't')
				&& (O[t + 2] == 'G' || O[t + 2] == 'g')
				&& ((O[t] == 'A' || O[t] == 'a')
					|| (O[t] == 'G' || O[t] == 'g')
					|| (O[t] == 'T' || O[t] == 't'))
			{
				alpha[S_STATE][t] = max_dbl;
				alpha[S_STATE][t + 1] = max_dbl;
				alpha[S_STATE][t + 2] = alpha[R_STATE][t - 1] - hmm.tr[TR_RS];
				path[S_STATE][t] = R_STATE as i8;
				path[S_STATE][t + 1] = S_STATE as i8;
				path[S_STATE][t + 2] = S_STATE as i8;

				let mut temp_alpha = alpha[E_STATE][t - 1] - hmm.tr[TR_ES];
				if temp_alpha < alpha[S_STATE][t + 2] {
					alpha[S_STATE][t + 2] = temp_alpha;
					path[S_STATE][t] = E_STATE as i8;
				}

				temp_alpha = alpha[E_STATE_1][t - 1] - hmm.tr[TR_ES1];
				if temp_alpha < alpha[S_STATE][t + 2] {
					alpha[S_STATE][t + 2] = temp_alpha;
					path[S_STATE][t] = E_STATE_1 as i8;
				}

				if O[t] == 'A' || O[t] == 'a' {
					alpha[S_STATE][t + 2] = alpha[S_STATE][t + 2] - log83;
				} else if O[t] == 'G' || O[t] == 'g' {
					alpha[S_STATE][t + 2] = alpha[S_STATE][t + 2] - (0.10 as f64).ln();
				} else if O[t] == 'T' || O[t] == 't' {
					alpha[S_STATE][t + 2] = alpha[S_STATE][t + 2] - log07;
				}

				/* adjustment based on probability distribution */
				let mut start_freq = 0.0;
				let mut sub_sum = 0.0;

				if t >= 30 {
					for i in 0..=60 {
						let i = i as isize - 30;
						if t as isize + i + 2 < len_seq as isize {
							start_freq -= hmm.tr_s[(i + 30) as usize][trinucleotide(
								&O[(t as isize + i) as usize],
								&O[(t as isize + i + 1) as usize],
								&O[(t as isize + i + 2) as usize],
							)];
						}
					}
				} else {
					let mut i = -1 * t as isize;
					while i <= 30 {
						if t as isize + i + 2 < len_seq as isize {
							sub_sum += hmm.tr_s[(i + 30) as usize][trinucleotide(
								&O[(t as isize + i) as usize],
								&O[(t as isize + i + 1) as usize],
								&O[(t as isize + i + 2) as usize],
							)];
						}
						i += 1;
					}
					sub_sum = sub_sum * 61.0 / (30.0 + t as f64 + 1.0);
					start_freq -= sub_sum;
				}

				let h_kd = hmm.s_dist[2]
					* (-1.0 * (start_freq - hmm.s_dist[1]).powi(2) / (2.0 * hmm.s_dist[0].powi(2)))
						.exp();
				let r_kd = hmm.s_dist[5]
					* (-1.0 * (start_freq - hmm.s_dist[4]).powi(2) / (2.0 * hmm.s_dist[3].powi(2)))
						.exp();
				let mut p_kd = h_kd / (h_kd + r_kd);
				if p_kd < 0.01 {
					p_kd = 0.01;
				} else if p_kd > 0.99 {
					p_kd = 0.99;
				}
				alpha[S_STATE][t + 2] = alpha[S_STATE][t + 2] - (p_kd).ln();
			}
		}

		/**********************************************/
		/* END' state                                 */
		/* originally start codon of genes in - strand */
		/**********************************************/
		if alpha[E_STATE_1][t] == 0.0 {
			alpha[E_STATE_1][t] = max_dbl;
			path[E_STATE_1][t] = NOSTATE;

			if t < len_seq - 2
				&& (O[t] == 'C' || O[t] == 'c')
				&& (O[t + 1] == 'A' || O[t + 1] == 'a')
				&& ((O[t + 2] == 'T' || O[t + 2] == 't')
					|| (O[t + 2] == 'C' || O[t + 2] == 'c')
					|| (O[t + 2] == 'A' || O[t + 2] == 'a'))
			{
				/* transition from frame6 */
				alpha[E_STATE_1][t + 2] = alpha[M6_STATE_1][t - 1] - hmm.tr[TR_GE];
				path[E_STATE_1][t] = M6_STATE_1 as i8;
				alpha[E_STATE_1][t] = max_dbl;
				alpha[E_STATE_1][t + 1] = max_dbl;
				path[E_STATE_1][t + 1] = E_STATE_1 as i8;
				path[E_STATE_1][t + 2] = E_STATE_1 as i8;

				if O[t + 2] == 'T' || O[t + 2] == 't' {
					alpha[E_STATE_1][t + 2] = alpha[E_STATE_1][t + 2] - log83;
				} else if O[t + 2] == 'C' || O[t + 2] == 'c' {
					alpha[E_STATE_1][t + 2] = alpha[E_STATE_1][t + 2] - (0.10 as f64).ln();
				} else if O[t + 2] == 'A' || O[t + 2] == 'a' {
					alpha[E_STATE_1][t + 2] = alpha[E_STATE_1][t + 2] - log07;
				}

				/* adjustment based on probability distribution */
				let mut start_freq = 0.0;

				let mut sub_sum = 0.0;

				if t >= 30 {
					for i in 0..=60 {
						let i = i as isize - 30;
						if t as isize + i + 2 < len_seq as isize {
							start_freq -= hmm.tr_e_1[(i + 30) as usize][trinucleotide(
								&O[(t as isize + i) as usize],
								&O[(t as isize + i + 1) as usize],
								&O[(t as isize + i + 2) as usize],
							)];
						}
					}
				} else {
					let mut i = -1 * t as isize;
					while i <= 30 {
						if t as isize + i + 2 < len_seq as isize {
							sub_sum += hmm.tr_e_1[(i + 30) as usize][trinucleotide(
								&O[(t as isize + i) as usize],
								&O[(t as isize + i + 1) as usize],
								&O[(t as isize + i + 2) as usize],
							)];
						}
						i += 1;
					}
					sub_sum = sub_sum * 61.0 / (30.0 + t as f64 + 1.0);
					start_freq -= sub_sum;
				}

				let h_kd = hmm.e1_dist[2]
					* (-1.0 * (start_freq - hmm.e1_dist[1]).powi(2)
						/ (2.0 * hmm.e1_dist[0].powi(2)))
					.exp();
				let r_kd = hmm.e1_dist[5]
					* (-1.0 * (start_freq - hmm.e1_dist[4]).powi(2)
						/ (2.0 * hmm.e1_dist[3].powi(2)))
					.exp();
				let mut p_kd = h_kd / (h_kd + r_kd);

				if p_kd < 0.01 {
					p_kd = 0.01;
				} else if p_kd > 0.99 {
					p_kd = 0.99;
				}
				alpha[E_STATE_1][t + 2] = alpha[E_STATE_1][t + 2] - (p_kd).ln();
			}
		}
		if num_N > 9 {
			for i in 0..NUM_STATE {
				if i != R_STATE {
					alpha[i][t] = max_dbl;
					path[i][t] = R_STATE as i8;
				}
			}
		}
	}

	/***********************************************************/
	/* backtrack array to find the optimal path                */
	/***********************************************************/

	let head_short: Vec<&str> = head.split_whitespace().collect();

	let mut out = String::with_capacity(5000);
	let mut aa = String::with_capacity(5000);
	let mut dna_output = String::with_capacity(5000);


	write!(out, ">{}\n", head_short[0]);

	/* find the state for O[N] with the highest probability */
	let mut prob = f64::INFINITY;
	for i in 0..NUM_STATE {
		if alpha[i][len_seq - 1] < prob {
			prob = alpha[i][len_seq - 1];
			vpath[len_seq - 1] = i;
		}
	}

	/* backtrack the optimal path */
	for t in (0..=(len_seq - 2)).rev() {
		vpath[t] = path[vpath[t + 1]][t + 1] as usize;
	}

	let mut codon_start = 0;
	let mut start_t: isize = -1;

	let mut codon = vec!['\0'; 4];
	let mut utr = vec!['\0'; 65];

	let mut dna = vec!['\0'; 300000];
	let mut dna1 = vec!['\0'; 300000];
	let mut dna_f = vec!['\0'; 300000];
	let mut dna_f1 = vec!['\0'; 300000];
	let mut protein = vec!['\0'; 100000];
	let mut insert = [0; 100];
	let mut delete = [0; 100];

	let mut start_orf = 0;
	let mut prev_match = 0;

	let mut insert_id = 0;
	let mut delete_id = 0;
	let mut end_t: usize;

	let mut final_score;
	let mut frame;

	let mut dna_start_t_withstop = 0;
	let mut dna_start_t = 0;

	for t in 0..len_seq {
		if codon_start == 0
			&& start_t < 0
			&& ((vpath[t] >= M1_STATE && vpath[t] <= M6_STATE)
				|| (vpath[t] >= M1_STATE_1 && vpath[t] <= M6_STATE_1)
				|| vpath[t] == S_STATE
				|| vpath[t] == S_STATE_1)
		{
			dna_start_t_withstop = t + 1;
			dna_start_t = t + 1;
			start_t = t as isize + 1;
			//introduce dna_start_t_withstop YY July 2018
		}

		if codon_start == 0
			&& (vpath[t] == M1_STATE
				|| vpath[t] == M4_STATE
				|| vpath[t] == M1_STATE_1
				|| vpath[t] == M4_STATE_1)
		{
			for i in 0..300000 {
				dna[i] = '\0';
				dna1[i] = '\0';
				dna_f[i] = '\0';
				dna_f1[i] = '\0';
			}
			for value in protein.iter_mut() {
				*value = '\0';
			}
			for value in insert.iter_mut(){
				*value = 0;
			}
			for value in delete.iter_mut(){
				*value = 0;
			}

			insert_id = 0;
			delete_id = 0;
			dna_id = 0;
			dna_f_id = 0;
			dna[dna_id] = O[t];
			dna_start_t_withstop = t + 1; //Ye April 21, 2016
			dna_start_t = t + 1;
			if vpath[t] == M1_STATE_1 || vpath[t] == M4_STATE_1 {
				if t > 2 {
					dna_start_t_withstop = t - 2;
				}
			}
			dna_f[dna_f_id] = O[t];
			//printf("Note start dna: t = %d, dna_id %d, dna_f_id %d, add %c\n", t, dna_id, dna_f_id, O[t]);
			start_orf = t + 1;
			prev_match = vpath[t];

			if vpath[t] < M6_STATE {
				codon_start = 1;
			} else {
				codon_start = -1;
			}
		} else if codon_start != 0
			&& (vpath[t] == E_STATE || vpath[t] == E_STATE_1 || t == len_seq - 1)
		{
			if vpath[t] == E_STATE || vpath[t] == E_STATE_1 {
				end_t = t + 3;
			} else {
				//end_t=t+1;
				/* FGS1.12 start: remove incomplete codon */
				let mut temp_t = t;
				while vpath[temp_t] != M1_STATE
					&& vpath[temp_t] != M4_STATE
					&& vpath[temp_t] != M1_STATE_1
					&& vpath[temp_t] != M4_STATE_1
				{
					dna_f[dna_f_id] = '\0';
					dna_f_id -= 1;

					dna[dna_id] = '\0';
					dna_id -= 1;

					temp_t -= 1;
				}
				end_t = temp_t; //??? YY July 2018
				      /* FGS1.12 end: remove incomplete codon */
			}

			if dna_id > gene_len {
				//these three lines moved here from outside of the loop above, YY July 23, 2018
				final_score = (alpha[vpath[(end_t - 4) as usize]][(end_t - 4) as usize]
					- alpha[vpath[(start_t + 2) as usize]][(start_t + 2) as usize])
					/ ((end_t as isize - start_t - 5) as f64);
				frame = start_orf % 3;
				if frame == 0 {
					frame = 3;
				}

				if codon_start == 1 {
					if start_t == dna_start_t as isize - 3 {
						//add complete start codon to dna, Ye April 21, 2016
						dna_start_t -= 3;
					}
					if refine {
						//add refinement of the start codons here, Ye, April 16, 2016
						let start_old = start_t;
						codon[0] = '\0';
						strncpy(&mut codon, O, (start_old - 1).try_into().unwrap(), 3);
						codon[3] = '\0';
						let mut s = 0;
						//find the optimal start codon within 30bp up- and downstream of start codon
						let mut e_save = 0.0;
						let mut s_save = 0; //initialization, YY July 25 2018

						let mut c = codon[0..3].iter().collect::<String>();
						while (!(c != "TAA"
							|| c != "TAG"
							|| c != "TGA"))
							&& (start_old - 1 - s - 35 >= 0)
						{
							if c != "ATG"
								|| c != "GTG"
								|| c != "TTG"
							{
								utr[0] = '\0';
								strncpy(
									&mut utr,
									O,
									(start_old - 1 - s - 30).try_into().unwrap(),
									63,
								);
								utr[63] = '\0';
								//printf("check s=%d, codon %s\n", s, codon);
								let mut freq_sum = 0.0;
								for j in 0..(strlen(&utr.iter().collect()) - 2) {
									let idx = trinucleotide(&utr[j], &utr[j + 1], &utr[j + 2]);
									freq_sum -= train.start[cg][j][idx];
									//printf("j=%d, key=%c%c%c %d, start %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->start[cg][j][idx]);
								}
								if s == 0 {
									e_save = freq_sum;
									s_save = 0;
								} else if freq_sum < e_save {
									e_save = freq_sum;
									s_save = -1 * s;
								} //positive chain, upstream s_save = -1 * 3
								 //printf("s=%d freq_sum %lf\n", s, freq_sum);
								 //getchar();
							}
							s += 3;
							codon[0] = '\0';
							strncpy(&mut codon, O, (start_old - 1 - s).try_into().unwrap(), 3);
							codon[3] = '\0';
							c = codon[0..3].iter().collect::<String>();
						}
						//update start_t YY July 2018
						if s_save != 0 {
							start_t = start_old + s_save;
							dna_start_t = (dna_start_t as isize + s_save) as usize;
						}
					}

					let dna_end_t = end_t;

					write!(
						out,
						"{}\t{}\t+\t{}\t{}\t",
						dna_start_t, dna_end_t, frame, final_score
					);
					write!(out, "I:");
					for i in 0..insert_id {
						write!(out, "{},", insert[i]);
					}
					write!(out, "\tD:");
					for i in 0..delete_id {
						write!(out, "{},", delete[i]);
					}
					write!(out, "\n");

					//update dna before calling get_protein, YY July 2018
					dna[0] = '\0';
					strncpy(
						&mut dna,
						O,
						dna_start_t - 1,
						dna_end_t - dna_start_t + 1,
					);
					dna[dna_end_t - dna_start_t + 1] = '\0';
					//end of update dna

					get_protein(&mut protein, &dna, true, wholegenome);
					/*
					if(!(strlen(protein) * 3 == strlen(dna) || strlen(protein) * 3 == strlen(dna) - 3)) {
					  printf("inconsistent protein/dna length: %d %d\n", strlen(protein), strlen(dna));
					  exit(0);
					}
					*/
					write!(
						aa,
						">{}_{}_{}_+\n",
						head_short[0], dna_start_t, dna_end_t
					);
					write!(
						dna_output,
						">{}_{}_{}_+\n",
						head_short[0], dna_start_t, dna_end_t
					);
					write!(aa, "{}\n", protein[0..veclen_str(&protein)].iter().collect::<String>());

					if !format {
						write!(dna_output, "{}\n", dna[0..(dna_end_t - dna_start_t + 1)].iter().collect::<String>());
					} else {
						write!(dna_output, "{}\n", dna_f[0..(dna_end_t - dna_start_t + 1)].iter().collect::<String>());
					}
				} else if codon_start == -1 {
					if refine {
						//add refinement of the start codons here, Ye, April 16, 2016
						let end_old = end_t; //reverse
						codon[0] = '\0';
						strncpy(&mut codon, O, (end_t - 3) as usize, 3);
						codon[3] = '\0';
						let mut s = 0;
						//find the optimal start codon within 30bp up- and downstream of start codon
						let mut e_save = 0.0;
						let mut s_save = 0; //initialization, YY July 25, 2018
						let mut c = codon[0..3].iter().collect::<String>();

						while (!(c != "TTA"
							|| c != "CTA"
							|| c != "TCA"))
							&& (end_old - 2 + s + 35 < len_seq)
						{
							if c != "CAT"
								|| c != "CAC"
								|| c != "CAA"
							{
								utr[0] = '\0';
								strncpy(&mut utr, O, (end_old - 3 + s - 30) as usize, 63);
								utr[63] = '\0';
								//printf("check s=%d, codon %s\n", s, codon);
								let mut freq_sum = 0.0;
								for j in 0..strlen(&utr.iter().collect()) - 2 {
									let idx = trinucleotide(&utr[j], &utr[j + 1], &utr[j + 2]);
									freq_sum -= train.stop1[cg][j][idx]; //stop1?? Ye, April 18, 2016
									 //printf("j=%d, key=%c%c%c %d, stop1 %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->stop1[cg][j][idx]);
								}
								if s == 0 {
									e_save = freq_sum;
									s_save = s;
								} else if freq_sum < e_save {
									e_save = freq_sum;
									s_save = s;
								} //negative chain, s_save = s, add, YY July 2018
							}
							s += 3;
							codon[0] = '\0';
							strncpy(&mut codon, O, (end_old - 3 + s) as usize, 3);
							codon[3] = '\0';
							c = codon[0..3].iter().collect::<String>();
						}
						//update end_t
						end_t = end_old + s_save;
					}

					let dna_end_t = end_t;
					write!(
						out,
						"{}\t{}\t-\t{}\t{}\t",
						dna_start_t_withstop, dna_end_t, frame, final_score
					);
					write!(out, "I:");
					for i in 0..insert_id {
						write!(out, "{},", insert[i]);
					}
					write!(out, "\tD:");
					for i in 0..delete_id {
						write!(out, "{},", delete[i]);
					}
					write!(out, "\n");

					//update dna before calling get_protein, YY July 2018
					//use dna_end_t & dna_start_w_withstop to avoid incomplete codons & include start/stop codons
					dna[0] = '\0';
					strncpy(
						&mut dna,
						O,
						dna_start_t_withstop - 1,
						dna_end_t - dna_start_t_withstop + 1,
					);
					dna[dna_end_t - dna_start_t_withstop + 1] = '\0';
					//end of update dna

					get_protein(&mut protein, &dna, false, wholegenome); //YY July 18, 2018, introduce adjust

					write!(
						aa,
						">{}_{}_{}_-\n",
						head_short[0], dna_start_t_withstop, dna_end_t
					);
					write!(
						dna_output,
						">{}_{}_{}_-\n",
						head_short[0], dna_start_t_withstop, dna_end_t
					);

					let dna1_out = get_rc_dna(&dna[0..((dna_end_t - dna_start_t_withstop + 1))].to_vec());
					let dna_f1_out = get_rc_dna_indel(&dna_f[0..((dna_end_t - dna_start_t_withstop + 1))].to_vec());
					write!(aa, "{}\n", protein[0..veclen_str(&protein)].iter().collect::<String>());
					if !format {
						write!(dna_output, "{}\n", dna1_out.iter().collect::<String>());
					} else {
						write!(dna_output, "{}\n", dna_f1_out.iter().collect::<String>());
					}
				}
			}
			codon_start = 0;
			start_t = -1;
			dna_id = 0;
			dna_f_id = 0;
		} else if codon_start != 0
			&& ((vpath[t] >= M1_STATE && vpath[t] <= M6_STATE)
				|| (vpath[t] >= M1_STATE_1 && vpath[t] <= M6_STATE_1))
			&& (vpath[t] as i8) - (prev_match as i8) < 6
		{
			let out_nt;
			if vpath[t] < prev_match {
				out_nt = vpath[t] + 6 - prev_match;
			} else {
				out_nt = vpath[t] - prev_match;
			}
			for kk in 0..out_nt {
				/* for deleted nt in reads */
				dna_id += 1;
				dna[dna_id] = 'N';
				//printf("dna_id %d, dna-len %d\n", dna_id, strlen(dna));
				dna_f_id += 1;
				dna_f[dna_f_id] = 'x';
				if kk > 0 {
					delete[delete_id] = (t + 1) as i8;
					delete_id += 1;
				}
			}
			dna[dna_id] = O[t];
			//printf("dna_id %d, add %d %c dna-len %d\n", dna_id, t, O[t], strlen(dna));
			dna_f[dna_f_id] = O[t];
			prev_match = vpath[t];
		} else if codon_start != 0
			&& ((vpath[t] >= I1_STATE && vpath[t] <= I6_STATE)
				|| (vpath[t] >= I1_STATE_1 && vpath[t] <= I6_STATE_1))
		{
			dna_f_id += 1;
			dna_f[dna_f_id] = O[t].to_lowercase().collect::<Vec<char>>()[0]; // TODO:  not so elegant
			insert[insert_id] = (t + 1) as i8;
			insert_id += 1;
		} else if codon_start != 0 && vpath[t] == R_STATE {
			/* for long NNNNNNNNN, pretend R state */
			codon_start = 0;
			start_t = -1;
			dna_id = 0;
			dna_f_id = 0;
		}
	};
	Output {
		out: out,
		aa: aa,
		dna: dna_output
	}
}

fn strlen(s: &String) -> usize {
	let mut result = 0;
	for c in s.chars() {
		if c != '\0' {
			result += 1;
		}
		else{
			break
		}
	}
	result
}


fn veclen_str(s: &Vec<char>) -> usize { // TODO: cleanup
	let mut result = 0;
	for c in s.iter() {
		if *c != '\0' {
			result += 1;
    }
    else{
      break;
    }
	}
	result
}

fn strncpy(target: &mut Vec<char>, source: &Vec<char>, start_pos: usize, len: usize) {
	let sourcelen = source.len();
	for i in start_pos..(start_pos + len) {
		if i < sourcelen {
			target[i - start_pos] = source[i];
		} else {
			target[i - start_pos] = '\0';
		}
	}
}

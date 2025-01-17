#pragma warning(disable:4996) //memcpy
#include "Needleman_Wunsch.h"

Needleman_Wunsch::Needleman_Wunsch(char* a, char* b, int gap_score, int match_score, int mismatch_score)
{
	this->gap_score = gap_score;
	this->match_score = match_score;
	this->mismatch_score = mismatch_score;

	this->aa = a;
	this->bb = b;

	int lena = strlen(a) + 1; // Help length variable ScoringMatrix lena - alignment
	int lenb = strlen(b) + 1; // Help length variable ScoringMatrix lenb - reference genome

	// Scoring Matrix Inicialization
	this->ScoringMatrix = new int *[lena];
	for (int i = 0; i < lena; i++) {
		this->ScoringMatrix[i] = new int[lenb];
	}
	
	// Fill first row and first column with default score
	this->ScoringMatrix[0][0] = 0;
	for (int i = 1; i < lena; i++) {
		ScoringMatrix[i][0] = i * this->gap_score;
	}

	for (int j = 1; j < lenb; j++) {
		ScoringMatrix[0][j] = j * this->gap_score;
	}

	/* 
	Fill Matrix
	*/
	for (int i = 1; i < lena; i++) {
		for (int j = 1; j < lenb; j++) {
			int score = 0;
			score = CalculateScore(i, j);
			this->matrix_max = score;
			this->i_max = i;
			this->j_max = j;
			this->ScoringMatrix[i][j] = score;
		}
	}

	/*
		Traceback
		diagonal: match/mismatch
		up:       gap in sequence 1
		left:     gap in sequence 2
	*/
	int m = this->i_max;
	int n = this->j_max;
	bool first = false;	//Help to indicate fist position and don't allow next iteration to rewrite it.
	while (m > 0 && n > 0)
	{
		int score = 0;
		if ((m == 1 || n == 1) && first == false) {
			first = true;
			this->i_min = m;
			this->j_min = n;
		}
		if (this->aa[m - 1] == this->bb[n - 1]) {
			score = this->match_score;
		}
		else {
			score = this->mismatch_score;
		}
		//Move diag
		if (m > 0 && n > 0 && ScoringMatrix[m][n] == ScoringMatrix[m - 1][n - 1] + score)
		{
			if ((this->aa[m - 1]) != (this->bb[n - 1])) {
				this->mismatch++;
			}
			this->cigar_str = "M" + this->cigar_str;
			m--; n--;
		}
		//Move left
		else if (n > 0 && ScoringMatrix[m][n] == ScoringMatrix[m][n - 1] + this->gap_score)
		{
			this->cigar_str = "D" + this->cigar_str;
			n--;
		}
		//Move up
		else //if (m > 0 && ScoringMatrix[m][n] == ScoringMatrix[m - 1][n] + this->gap_score)
		{
			this->cigar_str = "I" + this->cigar_str;
			m--;
		}
	}

	//CIGAR string computing
	this->cigar = strdup(this->cigar_str.c_str());
	CIGAR *cigar_cls = new CIGAR(this->cigar);
	this->cigar_final = cigar_cls->get_CIGAR();
	this->cigar_length = cigar_cls->get_CIGAR_length();
	delete cigar_cls;
}

int Needleman_Wunsch::CalculateScore(int i, int j)
{
	int max = 0;
	int diag_score = 0;
	int left_score = 0;
	int up_score = 0;
	//Set diag/left/up score
	if (this->aa[i - 1] == this->bb[j - 1]) {
		diag_score = this->ScoringMatrix[i - 1][j - 1] + this->match_score;
	}
	else {
		diag_score = this->ScoringMatrix[i - 1][j - 1] + this->mismatch_score;
	}
	left_score = this->ScoringMatrix[i][j - 1] + this->gap_score;
	up_score = this->ScoringMatrix[i-1][j] + this->gap_score;
	//Set maximum score from diag/left/up score
	if (diag_score > max) {
		max = diag_score;
	}
	if (left_score > max) {
		max = left_score;
	}
	if (up_score > max) {
		max = up_score;
	}
	return max;
}

int Needleman_Wunsch::get_first_pos() {
	return this->j_min;
}

int Needleman_Wunsch::get_last_pos() {
	return this->j_max;
}

string Needleman_Wunsch::get_cigar() {
	return this->cigar_final;
}

int Needleman_Wunsch::get_cigar_length() {
	return this->cigar_length;
}

int Needleman_Wunsch::get_matrix_max_score() {
	return this->matrix_max;
}

int Needleman_Wunsch::get_mismatch() {
	return this->mismatch;
}

Needleman_Wunsch::~Needleman_Wunsch()
{
	int lena = strlen(aa) + 1;
	int lenb = strlen(bb) + 1;
	for (int i = 0; i < lena; i++) {
		for (int j = 0; j < lenb; j++) {
		}
		delete[] ScoringMatrix[i];
	}
	ScoringMatrix = NULL;
	aa = NULL;
	bb = NULL;
	delete[] cigar;
	cigar = NULL;
}
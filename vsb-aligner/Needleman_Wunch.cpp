#pragma warning(disable:4996) //memcpy
#include "Needleman_Wunch.h"

Needleman_Wunch::Needleman_Wunch(char* a, char* b, int gap_score, int match_score, int mismatch_score)
{
	this->gap_score = gap_score;
	this->match_score = match_score;
	this->mismatch_score = mismatch_score;

	this->aa = a;
	this->bb = b;

	int lena = strlen(a) + 1; // Pomocná promenná pro ScoringMatrix lena - read
	int lenb = strlen(b) + 1; // Pomocná promenná pro ScoringMatrix lenb - reference genome

	// Scoring Matrix Inicialization
	this->ScoringMatrix = new int *[lena];
	for (int i = 0; i < lena; i++) {
		this->ScoringMatrix[i] = new int[lenb];
	}
	
	// Fill first row and first column with default score
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
			if (score > this->matrix_max)
			{
				this->matrix_max = score;
				this->i_max = i;
				this->j_max = j;
			}
			this->ScoringMatrix[i][j] = score;
		}
	}

	/*
		Traceback
		diagonal: match/mismatch
		up:       gap in sequence 1
		left:     gap in sequence 2
	*/
	stack<char> SCigar;
	int m = this->i_max;
	int n = this->j_max;
	while (m > 0 && n > 0)
	{
		int score = 0;
		if (aa[m - 1] == bb[n - 1]) {
			score = this->match_score;
		}
		else {
			score = this->mismatch_score;
		}
		if (m > 0 && n > 0 && ScoringMatrix[m][n] == ScoringMatrix[m - 1][n - 1] + score)
		{
			if ((this->aa[m - 1]) != (this->bb[n - 1])) {
				this->mismatch++;
			}
			SCigar.push('M');
			m--; n--;
		}
		else if (n > 0 && ScoringMatrix[m][n] > ScoringMatrix[m][n - 1] - this->gap_score)
		{
			SCigar.push('D');
			n--;
		}
		else
		{
			SCigar.push('I');
			m--;
		}
	}

	//reverse CIGAR string
	while (!SCigar.empty())
	{
		this->cigar_str += SCigar.top();
		SCigar.pop();
	}

	//CIGAR string computing
	this->cigar = strdup(this->cigar_str.c_str());
	CIGAR *cigar_cls = new CIGAR(this->cigar);
	this->cigar_final = cigar_cls->get_CIGAR();
	this->cigar_length = cigar_cls->get_CIGAR_length();
	delete cigar_cls;
}

int Needleman_Wunch::CalculateScore(int i, int j)
{
	int max = 0;
	int similarity;
	if (this->aa[i - 1] == this->bb[j - 1]) {
		similarity = this->match_score;
	}
	else {
		similarity = this->mismatch_score;
	}


	if ((ScoringMatrix[i - 1][j] + this->gap_score) >= (ScoringMatrix[i][j - 1] + this->gap_score)) {
		max = (ScoringMatrix[i - 1][j] + this->gap_score);
	}
	else {
		max = (ScoringMatrix[i][j - 1] + this->gap_score);
	}
	if ((ScoringMatrix[i - 1][j - 1] + similarity) >= max) {
		return (ScoringMatrix[i - 1][j - 1] + similarity);
	}
	else {
		return max;
	}
}

int Needleman_Wunch::get_first_pos() {
	return this->j_min;
}

int Needleman_Wunch::get_last_pos() {
	return this->j_max;
}

string Needleman_Wunch::get_cigar() {
	return this->cigar_final;
}

int Needleman_Wunch::get_cigar_length() {
	return this->cigar_length;
}

int Needleman_Wunch::get_matrix_max_score() {
	return this->matrix_max;
}

int Needleman_Wunch::get_mismatch() {
	return this->mismatch;
}


Needleman_Wunch::~Needleman_Wunch()
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
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
	string RetA, RetB;
	stack<char> SA, SB, SCigar;

	//char *xA, *xB, *xCigar;
	int lenc = strlen(a);
	int lend = strlen(b);
	if ((this->aa[lenc]) != (this->bb[lend])) {
		this->mismatch++;
	}
	SA.push(this->aa[lenc]);
	SB.push(this->bb[lend]);
	SCigar.push('M');
	while (lenc != 0 && lend != 0) // Change (lenc != 0 || lend != 0) to (lenc != 0 && lend != 0) to to filter out prefix with Inserts and Deletes.
	{
		if (lenc == 0)
		{
			SA.push('-');
			SB.push(bb[lend - 1]);
			SCigar.push('D');
			lend--;
		}
		else if (lend == 0)
		{
			SA.push(aa[lenc - 1]);
			SB.push('-');
			SCigar.push('I');
			lenc--;
		}
		else
		{
			int score = 0;
			if (aa[lenc - 1] == bb[lend - 1]) {
				score = this->match_score;
			}
			else {
				score = this->mismatch_score;
			}
			if (ScoringMatrix[lenc][lend] == ScoringMatrix[lenc - 1][lend - 1] + score)
			{
				if ((this->aa[lenc - 1]) != (this->bb[lend - 1])) {
					this->mismatch++;
				}
				SA.push(aa[lenc - 1]);
				SB.push(bb[lend - 1]);
				SCigar.push('M');
				lenc--; lend--;
			}
			else if (ScoringMatrix[lenc - 1][lend] > ScoringMatrix[lenc][lend - 1])
			{
				SA.push(aa[lenc - 1]);
				SB.push('-');
				SCigar.push('I');
				lenc--;
			}
			else
			{
				SA.push('-');
				SB.push(bb[lend - 1]);
				SCigar.push('D');
				lend--;
			}
		}
	}

	while (!SA.empty())
	{
		RetA += SA.top();
		RetB += SB.top();
		this->cigar_str += SCigar.top();
		SA.pop();
		SB.pop();
		SCigar.pop();
	}

	//CIGAR String computing
	this->cigar = new char[this->cigar_str.length()];
	strcpy(this->cigar, this->cigar_str.c_str());
	this->cigar = Cigar(this->cigar);
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

char* Needleman_Wunch::Cigar(char *a) {
	char* out;
	string test = "";
	int M = 0;
	int I = 0;
	int D = 0;
	for (int i = 0; i < strlen(a); i++) {
		if (a[i] == 'M') {
			M += 1;
			this->cigar_length += 1;
			if (I > 0) {
				test += to_string(I);
				test += 'I';
				I = 0;
			}
			else if (D > 0) {
				test += to_string(D);
				test += 'D';
				D = 0;
			}
		}
		else if (a[i] == 'I') {
			I += 1;
			this->cigar_length += 1;
			if (M > 0) {
				test += to_string(M);
				test += 'M';
				M = 0;
			}
			else if (D > 0) {
				test += to_string(D);
				test += 'D';
				D = 0;
			}
		}
		else if (a[i] == 'D') {
			D += 1;
			this->cigar_length += 1;
			if (M > 0) {
				test += to_string(M);
				test += 'M';
				M = 0;
			}
			else if (I > 0) {
				test += to_string(I);
				test += 'I';
				I = 0;
			}
		}
	}

	if (M > 0) {
		test += to_string(M);
		test += 'M';
		M = 0;
	}
	else if (I > 0) {
		test += to_string(I);
		test += 'I';
		I = 0;
	}
	else if (D > 0) {
		test += to_string(D);
		test += 'D';
		D = 0;
	}

	out = new char[test.length()];
	strcpy(out, test.c_str());

	return out;
}

int Needleman_Wunch::get_first_pos() {
	return this->j_min;
}

int Needleman_Wunch::get_last_pos() {
	return this->j_max;
}

char* Needleman_Wunch::get_cigar() {
	return this->cigar;
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
	this->gap_score = 0;
	this->match_score = 0;
	this->mismatch_score = 0;
	this->ScoringMatrix = 0;
	this->matrix_max = 0;
	this->i_max = 0, j_max = 0;
}
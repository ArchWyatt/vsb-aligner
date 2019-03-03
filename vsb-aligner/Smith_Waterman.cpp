#pragma warning(disable:4996) //memcpy
#include "Smith_Waterman.h"

Smith_Waterman::Smith_Waterman(char* a, char* b, int gap_score, int match_score, int mismatch_score)
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

	// Fill first row and first column with zero
	this->ScoringMatrix[0][0] = 0;
	for (int i = 1; i < lena; i++) {
		this->ScoringMatrix[i][0] = 0;
	}

	for (int j = 1; j < lenb; j++) {
		this->ScoringMatrix[0][j] = 0;
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

	/*	Traceback
	diagonal: match/mismatch
	up:       gap in sequence 1
	left:     gap in sequence 2
	*/
	stack<char> SCigar;
	int pos_i = this->i_max;
	int pos_j = this->j_max;
	int move = NextMove(pos_i, pos_j);
	while (move != 0) {
		if (move == 1) { // Diagonal move
			if ((this->aa[pos_i - 1]) != (this->bb[pos_j - 1])) {
				this->mismatch++;
			}
			SCigar.push('M');
			pos_i = pos_i - 1;
			pos_j = pos_j - 1;
		}
		else if (move == 2) { // UP move
			SCigar.push('I');
			pos_i = pos_i - 1;
		}
		else if (move == 3) { // LEFT move
			SCigar.push('D');
			pos_j = pos_j - 1;
		}
		else { // Error, program should not move 
			cout << "Chyba programu" << endl;
		}
		move = NextMove(pos_i, pos_j);
	}
	this->i_min = pos_i;
	this->j_min = pos_j;
	SCigar.push('M');

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

int Smith_Waterman::CalculateScore(int i, int j)
{
	int max = 0;
	int similarity;
	if (this->aa[i - 1] == this->bb[j - 1]) {
		similarity = this->match_score;
	}
	else {
		similarity = this->mismatch_score;
	}

	int diag_score = this->ScoringMatrix[i - 1][j - 1] + similarity;
	int up_score = this->ScoringMatrix[i - 1][j] + this->gap_score;
	int left_score = this->ScoringMatrix[i][j - 1] + this->gap_score;

	if (diag_score > max) {
		max = diag_score;
	}
	if (up_score > max) {
		max = up_score;
	}
	if (left_score > max) {
		max = left_score;
	}
	return max;
}

int Smith_Waterman::NextMove(int pos_i, int pos_j) {
	int diag = this->ScoringMatrix[pos_i - 1][pos_j - 1];
	int up = this->ScoringMatrix[pos_i - 1][pos_j];
	int left = this->ScoringMatrix[pos_i][pos_j - 1];
	if (diag >= up && diag >= left) { // Move diag
		if (diag != 0) { // 1 signals a DIAG move. 0 signals the end.
			return (int)1;
		}
		else {
			return 0;
		}
	}
	else if (up > diag && up >= left) {
		if (up != 0) { // 2 signals a UP move. 0 signals the end.
			return (int)2;
		}
		else {
			return 0;
		}
	}
	else if (left > diag && left > up) {
		if (left != 0) { // 3 signals a LEFT move. 0 signals the end.
			return (int)3;
		}
		else {
			return 0;
		}
	}
	else {
		cout << "Chyba programu" << endl;
		return 4;
	}
}

int Smith_Waterman::get_first_pos() {
	return this->j_min;
}

int Smith_Waterman::get_last_pos() {
	return this->j_max;
}

string Smith_Waterman::get_cigar() {
	return this->cigar_final;
}

int Smith_Waterman::get_cigar_length() {
	return this->cigar_length;
}

int Smith_Waterman::get_matrix_max_score() {
	return this->matrix_max;
}

int Smith_Waterman::get_mismatch() {
	return this->mismatch;
}

Smith_Waterman::~Smith_Waterman() {
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
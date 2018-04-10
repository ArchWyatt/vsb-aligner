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

	/* Fill Matrix

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


	//Print out of the matrix and fined max score in the matrix
	cout << "Testing String  : " << this->aa << endl;
	cout << "Reference String: " << this->bb << endl;
	cout << "Matrix:" << endl;
	for (int i = 0; i < lena; i++) {
		for (int j = 0; j < lenb; j++) {
			cout << right << setw(4) << ScoringMatrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "Max score in the matrix is " << this->matrix_max << endl;
	cout << "Matrix max score coordinates is i: " << this->i_max << " j: " << this->j_max << endl;

	/*	Traceback
	diagonal: match/mismatch
	up:       gap in sequence 1
	left:     gap in sequence 2
	*/
	string RetA, RetB;
	stack<char> SA, SB, SCigar;
	int pos_i = this->i_max;
	int pos_j = this->j_max;
	int move = NextMove(pos_i, pos_j);
	while (move != 0) {
		if (move == 1) { // Diagonal move
			SA.push(this->aa[pos_i - 1]);
			SB.push(this->bb[pos_j - 1]);
			SCigar.push('M');
			pos_i = pos_i - 1;
			pos_j = pos_j - 1;
		}
		else if (move == 2) { // UP move
			SA.push(this->aa[pos_i - 1]);
			SB.push('-');
			SCigar.push('I');
			pos_i = pos_i - 1;
		}
		else if (move == 3) { // LEFT move
			SA.push('-');
			SB.push(this->bb[pos_j - 1]);
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
	SA.push(this->aa[pos_i - 1]);
	SB.push(this->bb[pos_j - 1]);
	SCigar.push('M');

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

	cout << "AlignmentA          : " << RetA << endl;
	cout << "AlignmentB Reference: " << RetB << endl;
	cout << "Cigar               : " << this->cigar << endl;

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

char* Smith_Waterman::Cigar(char *a) {
	char* out;
	string test = "";
	int M = 0;
	int I = 0;
	int D = 0;
	for (int i = 0; i < strlen(a); i++) {
		if (a[i] == 'M') {
			M += 1;
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

int Smith_Waterman::get_first_pos() {
	return this->j_min;
}

int Smith_Waterman::get_last_pos() {
	return this->j_max;
}

char* Smith_Waterman::get_cigar() {
	return this->cigar;
}

int Smith_Waterman::get_matrix_max_score() {
	return this->matrix_max;
}


Smith_Waterman::~Smith_Waterman()
{
	this->gap_score = 0;
	this->match_score = 0;
	this->mismatch_score = 0;
	this->ScoringMatrix = 0;
	this->matrix_max = 0;
	this->i_max = 0, j_max = 0;
}
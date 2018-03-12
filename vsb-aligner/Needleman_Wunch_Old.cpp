#include "Needleman_Wunch_Old.h"

Needleman_Wunch_Old::Needleman_Wunch_Old(string a, string b, int gap_score, int match_score, int mismatch_score)
{
	this->a = a;
	this->b = b;
	this->gap_score = gap_score;
	this->match_score = match_score;
	this->mismatch_score = mismatch_score;

	int lena = this->a.length() + 1;
	int lenb = this->b.length() + 1;

	// Scoring Matrix Inicialization

	int **ScoringMatrix;
	ScoringMatrix = new int *[lena];
	for (int i = 0; i < lena; i++) {
		ScoringMatrix[i] = new int[lenb];
	}
	ScoringMatrix[0][0] = 0;

	// Fill first row and first column with default score
	for (int i = 1; i < lena; i++) {
		ScoringMatrix[i][0] = i * this->gap_score;
	}

	for (int j = 1; j < lenb; j++) {
		ScoringMatrix[0][j] = j * this->gap_score;
	}


	// Fill Matrix
	for (int i = 1; i < lena; i++) {
		for (int j = 1; j < lenb; j++) {
			//int score = (a[i - 1] == b[j - 1]) ? match_score : mismatch_score;
			int score = 0;
			if (a[i - 1] == b[j - 1]) {
				score = this->match_score;
			}
			else {
				score = this->mismatch_score;
			}
			//Max function replacement.
			/*
			int helper = 0;
			if ((ScoringMatrix[i - 1][j] + this->gap_score) >= (ScoringMatrix[i][j - 1] + this->gap_score)) {
				helper = (ScoringMatrix[i - 1][j] + this->gap_score);
			}
			else {
				helper = (ScoringMatrix[i][j - 1] + this->gap_score);
			}
			if ((ScoringMatrix[i - 1][j - 1] + score) >= helper) {
				ScoringMatrix[i][j] = (ScoringMatrix[i - 1][j - 1] + score);
			}
			else {
				ScoringMatrix[i][j] = helper;
			}
			*/
			ScoringMatrix[i][j] = max(ScoringMatrix[i - 1][j - 1] + score, max(ScoringMatrix[i - 1][j] + this->gap_score, ScoringMatrix[i][j - 1] + this->gap_score));
		}
	}


	//Print out of the matrix
	cout << "Testing String  : " << a << endl;
	cout << "Reference String: " << b << endl;
	cout << "Matrix:" << endl;
	for (int i = 0; i < lena; i++) {
		for (int j = 0; j < lenb; j++) {
			cout << right << setw(4) << ScoringMatrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;


	//Traceback
	string RetA, RetB, RetCigar;
	stack<char> SA, SB, CIGAR;
	int lenc = this->a.length();
	int lend = this->b.length();

	while (lenc != 0 || lend != 0)
	{
		if (lenc == 0)
		{
			SA.push('-');
			SB.push(b[lend - 1]);
			lend--;
		}
		else if (lend == 0)
		{
			SA.push(a[lenc - 1]);
			SB.push('-');
			lenc--;
		}
		else
		{
			int score = 0;
			if (a[lenc - 1] == b[lend - 1]) {
				score = this->match_score;
			}
			else {
				score = this->mismatch_score;
			}
			if (ScoringMatrix[lenc][lend] == ScoringMatrix[lenc - 1][lend - 1] + score)
			{
				SA.push(a[lenc - 1]);
				SB.push(b[lend - 1]);
				lenc--; lend--;
			}
			else if (ScoringMatrix[lenc - 1][lend] > ScoringMatrix[lenc][lend - 1])
			{
				SA.push(a[lenc - 1]);
				SB.push('-');
				lenc--;
			}
			else
			{
				SA.push('-');
				SB.push(b[lend - 1]);
				lend--;
			}
		}
	}

	while (!SA.empty())
	{
		RetA += SA.top();
		RetB += SB.top();
		SA.pop();
		SB.pop();
	}

	cout << "AlignmentA: " << RetA << endl;
	cout << "AlignmentB: " << RetB << endl;
}


Needleman_Wunch_Old::~Needleman_Wunch_Old()
{
	this->a = "";
	this->b = "";
	this->gap_score = 0;
	this->match_score = 0;
	this->mismatch_score = 0;
}
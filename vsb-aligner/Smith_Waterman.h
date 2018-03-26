#include <iostream>     // std::cout, std::endl
#include <string>
#include <algorithm>    // std::max
#include <iomanip>      // std::setw
#include <stack>

using namespace std;

class Smith_Waterman
{
private:
	string a = "";				//Sequence 1
	string b = "";				//Sequence 2
	char *aa;					//Sequence 1 in char array
	char *bb;					//Sequence 2 in char array
	int gap_score = 0;
	int match_score = 0;
	int mismatch_score = 0;
	int **ScoringMatrix;		//Score Matrix
	int matrix_max = 0;			//The highest score in matrix.
	int i_max = 0, j_max = 0;	//The row and columbn coordinates of the highest score in matrix.

public:
	Smith_Waterman(string a, string b, int gap_score, int match_score, int mismatch_score);
	int CalculateScore(int i, int j);
	int NextMove(int pos_i, int pos_j);
	~Smith_Waterman();
};
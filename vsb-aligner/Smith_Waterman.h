#include <iostream>     // std::cout, std::endl
#include <string>
#include <algorithm>    // std::max
#include <iomanip>      // std::setw
#include <stack>

using namespace std;

class Smith_Waterman
{
private:
	char *aa;					//Read
	char *bb;					//Reference genome
	char *cigar;				//Temporary cigar before final version
	string cigar_str = "";		//Cigar string computed in Smith_Waterman
	int gap_score = 0;
	int match_score = 0;
	int mismatch_score = 0;
	int **ScoringMatrix;		//Score Matrix
	int matrix_max = 0;			//The highest score in matrix.
	int i_min = 0, j_min = 0;	//The row and column coordinates of the lowest score in matrix. (j is the first position of aligned string)
	int i_max = 0, j_max = 0;	//The row and columbn coordinates of the highest score in matrix. (j is the last position of aligned string)

public:
	Smith_Waterman(char* a, char* b, int gap_score, int match_score, int mismatch_score);
	int CalculateScore(int i, int j);
	int NextMove(int pos_i, int pos_j);
	~Smith_Waterman();
	int get_first_pos();
	int get_last_pos();
	char* Cigar(char *a);
	char* get_cigar();
	int get_matrix_max_score();

};
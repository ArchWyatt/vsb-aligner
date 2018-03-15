#include <iostream>     // std::cout, std::endl
#include <string>
#include <algorithm>    // std::max
#include <iomanip>      // std::setw
#include <stack>

using namespace std;

class Needleman_Wunch_Old
{
private:
	string a = "";
	string b = "";
	int gap_score = 0;
	int match_score = 0;
	int mismatch_score = 0;

public:
	Needleman_Wunch_Old(string a, string b, int gap_score, int match_score, int mismatch_score);
	~Needleman_Wunch_Old();
};
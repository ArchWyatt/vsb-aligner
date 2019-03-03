#pragma once

#include <cstdio>
#include <string>

using namespace std;

class CIGAR
{
private:
	char *cigar;			//Cigar long version string
	int cigar_length = 0;	//Cigar Length
	

public:
	CIGAR() {};
	CIGAR(char *cigar) {
		this->cigar = cigar;
	}
	string get_CIGAR() {
		string temp_str = "";
		int M = 0;
		int I = 0;
		int D = 0;
		bool first_M = false;	//Start CIGAR string with M match/mismatch
		for (int i = 0; i < strlen(this->cigar); i++) {
			if (this->cigar[i] == 'M') {
				if (first_M == false) {
					first_M = true;
				}
				M += 1;
				this->cigar_length += 1;
				if (I > 0) {
					temp_str += to_string(I);
					temp_str += 'I';
					I = 0;
				}
				else if (D > 0) {
					temp_str += to_string(D);
					temp_str += 'D';
					D = 0;
				}
			}
			else if ((this->cigar[i] == 'I') && (first_M == true)) {
				I += 1;
				this->cigar_length += 1;
				if (M > 0) {
					temp_str += to_string(M);
					temp_str += 'M';
					M = 0;
				}
				else if (D > 0) {
					temp_str += to_string(D);
					temp_str += 'D';
					D = 0;
				}
			}
			else if ((this->cigar[i] == 'D') && (first_M == true)) {
				D += 1;
				this->cigar_length += 1;
				if (M > 0) {
					temp_str += to_string(M);
					temp_str += 'M';
					M = 0;
				}
				else if (I > 0) {
					temp_str += to_string(I);
					temp_str += 'I';
					I = 0;
				}
			}
		}

		if (M > 0) {
			temp_str += to_string(M);
			temp_str += 'M';
			M = 0;
		}
		/* We need only M match/mismatch at the end of the CIGAR string
		else if (I > 0) {
			temp_str += to_string(I);
			temp_str += 'I';
			I = 0;
		}
		else if (D > 0) {
			temp_str += to_string(D);
			temp_str += 'D';
			D = 0;
		}
		*/
		return temp_str;
	}

	int get_CIGAR_length() {
		return this->cigar_length;
	}

	~CIGAR() {};
};
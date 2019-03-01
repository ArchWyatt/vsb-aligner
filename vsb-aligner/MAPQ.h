#pragma once

#include <cstdio>

using namespace std;

class MAPQ
{
private:
	int k = 0;		//Mismatch
	char *quality;	//Quality array

public:
	MAPQ() {};
	MAPQ(int mismatch, char* quality) {
		this->k = mismatch;
		this->quality = quality;
	}
	int get_MAPQ() {
		int counter = 0;//Character counter
		int aq = 0;		//Average quality
		int temp = 0;	//Temporary MAPQ
		int str_len = strlen(this->quality);
		for (int i = 0; i < str_len; i++) {
			aq = aq + (int)quality[i] - 33;
			counter++;
		};
		aq = (aq / counter);
		temp = (4 + (3 - (this->k)) * (aq - 14));
		if (temp > 60) {
			return 60;
		}
		return temp;
	};
	~MAPQ() {};
};
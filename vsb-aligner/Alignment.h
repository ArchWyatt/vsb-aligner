#pragma once

#include "Definitions.h"

class Alignment
{
public:
	Alignment() {};
	Alignment(char* chrom, u_int c_pos) {
		chromosome = chrom;
		pos = c_pos;
	};

	~Alignment() {};

	/* c-string ukonceny nulou */
	char* chromosome;
	/* pozice zacatku */
	u_int pos;

	/* Sekce: Martin Kubala */
	char* cigar;
	u_int cigar_length;
	u_int score;
	/* Attribute for output, if it is true, alignment will be printed out */
	bool available = true;
	/* Mapping quality	*/
	u_int MAPQ = 0;
	/* FLAG	*/
	u_int FLAG = 0;
	/* Mark of alignment duplicity, to know, that there is same alignment */
	bool alternative = false;
};

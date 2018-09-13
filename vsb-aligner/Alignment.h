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
	/* proměnná pro output, v případě, že stejný aligment s vyšším score byl už zapsán, tak zde bude nastavena hodnota na false a už nebude dále porovnáván */
	bool available;
	/*
		Mapping quality
	*/
};

#pragma once

#include "Definitions.h"

class Alignment
{
public:
	Alignment() {};
	~Alignment() {};

	/* c-string ukonceny nulou */
	char* chromosome;
	/* pozice zacatku */
	u_int pos;

	/* Sekce: Martin Kubala*/
	/*
		Mapping quality
		CIGAR string
	*/
};

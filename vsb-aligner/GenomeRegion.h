#pragma once

#include "Definitions.h"

class GenomeRegion
{
public:
	GenomeRegion(char* chrom_id, u_int bases_no, u_int bases_s, u_int bases_in_line, u_int bytes_in_line) {
		chromosome_id = chrom_id;
		bases_number = bases_no;
		bases_start = bases_s;
		bases_per_line = bases_in_line;
		bytes_per_line = bytes_in_line;
		bases = NULL;
	};

	~GenomeRegion() {
		delete[] chromosome_id;

		if (bases != NULL) {
			delete[] bases;
		}
	};

	/*
		PROPERTIES
	*/
	char* chromosome_id;
	u_int bases_number;
	u_int bases_start;
	char* bases;
	u_int bases_per_line;
	u_int bytes_per_line;

	void SetBases(char* new_bases) {
		bases = new_bases;
	};
};
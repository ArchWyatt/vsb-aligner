#pragma once

/*
	Class representing one read
*/

#include "Alignment.h"
#include "Definitions.h"
#include "List.h"

static const u_char FORWARD_READ = 0;
static const u_char REVERSE_READ = 1;
static const u_char UNKNOWN_READ = 0xff;

class Read
{

public:
	Read() {		
	};

	Read(char* desc, char* seq, char* comm, char* qual, u_char type) {
		sequence = seq;
		descriptor = desc;		
		quality = qual;
		read_type = type;
		read_idx1 = 0;
		read_idx2 = 0;
		paired_read = NULL;	
	};

	~Read() {
		delete[] sequence;
		delete[] descriptor;
		delete[] name;
		delete[] quality;		
	};

	Read* paired_read;

	char* sequence;
	char* descriptor;
	char* name;	
	char* quality;
	u_char read_type;

	u_int read_idx1;
	u_int read_idx2;
		
	/*
		Pole alignmentu.
	*/
	List<Alignment*>* alignments;
	/*
		Delka pole alignmentu.
	*/
	u_char aligment_len;

	/*
		The expected length of the sequence based on the cigar and read beg
	*/
	u_int expected_len;
};
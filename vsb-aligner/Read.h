#pragma once

/*
	Class representing one read
*/

#include <fstream>

#include "Alignment.h"
#include "Definitions.h"
#include "List.h"

using namespace std;

static const u_char FORWARD_READ = 0;
static const u_char REVERSE_READ = 1;
static const u_char UNKNOWN_READ = 0xff;

class Read
{

public:
	Read() {		
		paired_read = NULL;
		alignments = new List<Alignment*>();
	};

	Read(char* desc, char* seq, char* comm, char* qual, u_char type) {
		sequence = seq;
		descriptor = desc;		
		quality = qual;
		read_type = type;
		read_idx1 = 0;
		read_idx2 = 0;

		paired_read = NULL;
		alignments = new List<Alignment*>();
	};

	~Read() {
		delete[] sequence;
		delete[] descriptor;
		delete[] name;
		delete[] quality;

		delete[] paired_read;
		delete alignments;
	};

	Read* paired_read;

	char* sequence;
	char* descriptor;
	char* name;	
	char* comment;
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

	/* to parse ids for pairing of reads */
	void ParseIDs() {
		u_int len = strlen(descriptor);
		u_int space_pos = 0;
		while (descriptor[space_pos] != ' ') space_pos++;

		name = new char[space_pos + 1];
		name[space_pos] = 0;
		memcpy(name, descriptor, sizeof(char)*space_pos);

		u_int first_id_pos = space_pos;
		while (descriptor[first_id_pos] != ':') first_id_pos--;
		first_id_pos += 1;

		u_int first_id_len = space_pos - first_id_pos + 1;
		char* first_id = new char[first_id_len + 1];
		first_id[first_id_len] = 0;
		memcpy(first_id, descriptor + first_id_pos, sizeof(char)*first_id_len);
		read_idx1 = atoi(first_id);
		delete[] first_id;

		u_int second_id_end = first_id_pos - 2;
		u_int second_id_sta = second_id_end;

		while (descriptor[second_id_sta] != ':') second_id_sta--;
		second_id_sta++;
		u_int second_id_len = second_id_end - second_id_sta + 1;
		char* second_id = new char[second_id_len + 1];
		second_id[second_id_len] = 0;
		memcpy(second_id, descriptor + second_id_sta, sizeof(char)*second_id_len);
		read_idx2 = atoi(second_id);
		delete[] second_id;
	};
};
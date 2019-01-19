#pragma once

/* type shortcuts */
typedef unsigned int u_int;
typedef unsigned char u_char;

#ifndef NULL
#define NULL 0
#endif

#ifndef ALIGNER_CONSTANTS
#define ALIGNER_CONSTANTS

#endif

#ifndef SA_CONSTANTS
#define SA_CONSTANTS

// SA = {$,A,C,G,N,T}
const u_int SA_ALFA_LEN = 6;

#endif

/*
Contains various parameters options for computation.
General syntax
- variable name => class_property = value
*/

class ProgOptions
{
public:
	ProgOptions() {
		/* Aligner setup - DEFAULTS*/
		/* Section: Martin Kubala */
		ID = "vsb";
		PN = "vsb-aligner";
		VN = "0.1";
		T = 0;

		//Smith Waterman score
		gap_score = -2;
		match_score = 3;
		mismatch_score = -3;

		range_prefix = 0;
		range_suffix = 0; 
	};

	~ProgOptions() {
	};	
	char* ID; // Program ID
	char* PN; // Program name
	char* VN; // program version
	int T; //Throw away readds with score less then specified number
	int gap_score;
	int match_score;
	int mismatch_score;
	u_int range_prefix;//Reference genom prefix n number of characters what will be included in calculation.
	u_int range_suffix;//Reference genom suffix n number of characters what will be included in calculation.
};

class ProgInfo
{
public:
	ProgInfo() {
		options = new ProgOptions();

		fq_F = NULL;
		fq_R = NULL;		
		genome_path = NULL;		
		sam_file = NULL;

		aligner_min_split_size = 20;
	};

	~ProgInfo() {		
		delete options;
	};

	/* Public properties */
	
	/* File with reads:
		- fq_F - forward read
		- fq_R - reverse read
	*/
	char* fq_F;
	char* fq_R;

	/* Manifest file */
	char* manifest_file;

	/* Input genome file */
	char* genome_path;
	
	/* SAM filepath*/
	char* sam_file;

	/* Computation related options */
	ProgOptions* options;

	/* Aligner */
	/* Pigeon hole minimum split size */
	u_char aligner_min_split_size;
};
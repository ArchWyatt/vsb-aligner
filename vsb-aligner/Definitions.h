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
	};

	~ProgOptions() {
	};	
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
};
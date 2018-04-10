#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>

#include "Definitions.h"
#include "Genome.h"
#include "List.h"
#include "Read.h"


using namespace std;

class Aligner
{
protected:
	ProgInfo* m_prog_info;

	List<Read>* forward_reads;
	List<Read>* reverse_reads;

	Genome* m_genome;
	
	void InitializeReads(char* reads_path, u_char read_type);

	void AlignOneRead(GenomeRegion* chromosome, SuffixArray* sa, Read* r);

public:
	Aligner();
	Aligner(ProgInfo* prog_info, Genome* genome);
	~Aligner();
	
	/* pair reads */
	void PairReads();

	/* each intervals read sequence without the first primer will be aligned */
	void AlignReads();

	List<Read>* Reads();

	static int Compare(const void* a, const void* b);
};
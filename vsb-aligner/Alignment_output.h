#pragma once

#include "Definitions.h"
#include "Alignment.h"

class Alignment_output
{
public:
	Alignment_output() {};
	Alignment_output(u_int idx, string namex, u_int FLAGx, char* RNAMEx, u_int POSx, u_int MAPQx, string CIGARx, char* RNEXTx, u_int PNEXTx, u_int TLENx, char* SEQx, char* QUALx, u_int scorex,bool available, bool alternative, bool top) {

		this->id = idx;
		this->NAME = namex;
		this->FLAG = FLAGx;
		this->RNAME = RNAMEx;
		this->POS = POSx;
		this->MAPQ = MAPQx;
		this->CIGAR = CIGARx;
		this->RNEXT = RNEXTx;
		this->PNEXT = PNEXTx;
		this->TLEN = TLENx;
		this->SEQ = SEQx;
		this->QUAL = QUALx;
		this->score = scorex;
		this->available = available;
		this->alternative = alternative;
		this->top = top;
	};

	~Alignment_output() {};

	/* Attributes */
	u_int id;
	string NAME;
	u_int FLAG;
	char* RNAME;
	u_int POS;
	u_int MAPQ;
	string CIGAR;
	char* RNEXT;
	u_int PNEXT;
	int TLEN;
	char* SEQ;
	char* QUAL;
	u_int score;
	bool available;
	bool alternative;
	bool top;
};
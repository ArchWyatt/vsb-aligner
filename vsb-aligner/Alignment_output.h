#pragma once

#include "Definitions.h"
#include "Alignment.h"

class Alignment_output
{
public:
	Alignment_output() {};
	Alignment_output(u_int idx, char* namex, u_int FLAGx, char* RNAMEx, u_int POSx, u_int MAPQx, char* CIGARx, char* RNEXTx, u_int PNEXTx, u_int TLENx, char* SEQx, char* QUALx, u_int scorex) {

		this->id = idx;
		this->QNAME = namex;
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
		this->available = true;



	};

	~Alignment_output() {};

	/* Attributes */
	u_int id;
	char* QNAME;
	u_int FLAG;
	char* RNAME;
	u_int POS;
	u_int MAPQ;
	char* CIGAR;
	char* RNEXT;
	u_int PNEXT;
	int TLEN;
	char* SEQ;
	char* QUAL;
	u_int score;
	bool available;

};
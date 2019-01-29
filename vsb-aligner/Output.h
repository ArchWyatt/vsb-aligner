#pragma once
#include <fstream>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>

#include "Alignment.h"
#include "Alignment_output.h"
#include "Definitions.h"
#include "Genome.h"
#include "GenomeRegion.h"
#include "List.h"
#include "Read.h"

using namespace std;

class Output {
private:
	ofstream ofs;
	char* sam_file = NULL;
	int score_MAX = 0;
	List<Alignment_output>* out_alignments = new List<Alignment_output>();
public:
	Output(char* output_file){
		this->sam_file = output_file;
		int file_check = remove(output_file); //If file exists, return code is 0. We have it here to be sure, that file will be deleted.
	}

	/*
	Print the head into the txt file
	*/
	void print_head(List<GenomeRegion>* gen){
		ListIterator<GenomeRegion> iterator(gen->First());
		this->ofs.open(this->sam_file, ios::out);
		while (iterator.Current() != NULL) {
			GenomeRegion* temp = iterator.Current()->Value();
			this->ofs << "@SQ SN:" << temp->chromosome_id << "\tLN:" << temp->bases_number << "\n";
			iterator.Next();
		}
		this->ofs.close();
	}

	/*
	Print program info into the txt file
	*/
	void print_program_info(char* ID, char* PN, char* VN, int T, char* r1, char* r2, char* genome) {
		this->ofs.open(this->sam_file, ios::out | ios::app);
		this->ofs << "@PG ID:" << ID << "\tPN:" << PN << "\tVN:" << VN << "\t-T " << T << "\t" << r1 << "\t" << r2 << "\t" << genome << "\n";
		this->ofs.close();
	}


	/*
	Preparing List of alignments output format for better handling 
	*/
	void output_prepare(List<Read>* reads, int T) {
		u_int id = 0;
		ListIterator<Read> iterator(reads->First());
		while (iterator.Current() != NULL) {
			Read* r = iterator.Current()->Value();
			ListIterator<Alignment> a_iterator(r->alignments->First());
			if (r->paired_read != NULL) {
				int temp = 0; // 
				Read* r2 = r->paired_read;
				ListIterator<Alignment> b_iterator(r2->alignments->First());
				while (a_iterator.Current() != NULL) {
					Alignment* a = a_iterator.Current()->Value();
					if ((a->score > T) && (a->available == true)) {
						u_int FLAG = 0;
						char* RNEXT;
						u_int PNEXT = 0;
						int TLEN = 0;
						/*
						RNEXT calculating part
						*/
						string x1 = r->name;
						string x2 = r->paired_read->name;
						if (x1 == x2) {
							RNEXT = "=";
						}
						else {
							RNEXT = "*";
						}
						/*
						PNEXT, TLEN calculating part
						*/
						if (b_iterator.Current() != NULL) {
							Alignment* b = b_iterator.Current()->Value();
							PNEXT = b->pos;
							if (b->pos >= a->pos) {
								temp = b->pos - a->pos;
								temp += b->cigar_length;
							}
							else {
								temp = a->pos - b->pos;
								temp += a->cigar_length;
							}
							TLEN = temp + 1;
							b_iterator.Next();
						}
						else {
							PNEXT = 0;
							TLEN = 0;
						}
						temp = 0;
						/*
						Store alignment into output list
						*/
						this->out_alignments->Append(new Alignment_output(id, r->name, FLAG, a->chromosome, a->pos, a->MAPQ, a->cigar, RNEXT, PNEXT, TLEN, r->sequence, r->quality, a->score));
						id++;
					}
					a_iterator.Next();
				}

				a_iterator = r->alignments->First();
				b_iterator = r2->alignments->First();
				while (b_iterator.Current() != NULL) {
					Alignment* b = b_iterator.Current()->Value();
					if ((b->score > T) && (b->available == true)) {
						u_int FLAG = 0;
						char* RNEXT;
						u_int PNEXT = 0;
						int TLEN = 0;
						/*
						RNEXT calculating part
						*/
						string x1 = r2->name;
						string x2 = r2->paired_read->name;
						if (x1 == x2) {
							RNEXT = "=";
						}
						else {
							RNEXT = "*";
						}
						/*
						PNEXT, TLEN calculating part
						*/
						if (a_iterator.Current() != NULL) {
							Alignment* a = a_iterator.Current()->Value();
							PNEXT = a->pos;
							if (a->pos >= b->pos) {
								temp = a->pos - b->pos;
								temp += a->cigar_length;

							}
							else {
								temp = b->pos - a->pos;
								temp += b->cigar_length;
							}
							TLEN = -temp-1;
							a_iterator.Next();
						}
						else {
							PNEXT = 0;
							TLEN = 0;
						}
						temp = 0;
						/*
						Store alignment into output list
						*/
						this->out_alignments->Append(new Alignment_output(id, r2->name, FLAG, b->chromosome, b->pos, b->MAPQ, b->cigar, RNEXT, PNEXT, TLEN, r2->sequence, r2->quality, b->score));
						id++;
					}
					b_iterator.Next();
				}

			}

			iterator.Next();
		}
	}

	/*
	Filter the alignments by the highes score of the alingments.
	*/
	void output_top_score_filtering() {
		//Zde bude protrizeni vyustupu kde bude mozne vypsat pouze alignmenty s nejvyssim score, pripadne alignmenty ktera budou mit spolecne nejvyssi score.
		//Vector bude pro alignmenty jenž jsou si rovny, bude se ukladat id prvku v listu a pokud bude prvek vyšší, tak se id použije k zneplatnění těchto alignmentu v samostatném cyklu dokud nebude vektor prázdný, takže while a následně bude porovnávaní nejvyššího alignmentu pokračovat.

	}
	
	/*
	Print alignments data into the file
	*/
	void print_output_data() {
		this->ofs.open(this->sam_file, ios::out | ios::app);
		ListIterator<Alignment_output> iterator(out_alignments->First());
		while (iterator.Current() != NULL) {
			Alignment_output* out = iterator.Current()->Value();
			if (out->available == true) {
				this->ofs << out->QNAME << "\t";
				this->ofs << out->FLAG << "\t";
				this->ofs << out->RNAME << "\t";
				this->ofs << out->POS << "\t";
				this->ofs << out->MAPQ << "\t";
				this->ofs << out->CIGAR << "\t";
				this->ofs << out->RNEXT << "\t";
				this->ofs << out->PNEXT << "\t";
				this->ofs << out->TLEN << "\t";
				this->ofs << out->SEQ << "\t";
				this->ofs << out->QUAL << "\t";
				this->ofs << "\n";
			}
			iterator.Next();
		}
		this->ofs.close();
	}

	void print_output_data_old(List<Read>* reads, int T) {
		//Tady se upravi pouze vypis z upraveneho listu
		this->ofs.open(this->sam_file, ios::out | ios::app);
		ListIterator<Read> iterator(reads->First());
		while (iterator.Current() != NULL) {
			Read* r = iterator.Current()->Value();
			ListIterator<Alignment> a_iterator(r->alignments->First());

			//Version only with paired read - for this work all of the reads is in the pair
			if (r->paired_read != NULL) {
				int temp = 0; // 
				Read* r2 = r->paired_read;
				ListIterator<Alignment> b_iterator(r2->alignments->First());

				//while (a_iterator.Current() != NULL && b_iterator.Current() != NULL) {
				while (a_iterator.Current() != NULL) {
					Alignment* a = a_iterator.Current()->Value();
					//Alignment* b = b_iterator.Current()->Value();
					/*
					//Start of the SCORE control to find best match or equal matches.
					if (a->available == true) {
						this->score_MAX = a->score;

						Here is the compare with all alignments, rules below
						- Compare this alignment with all Read & paired read alignments.
						- Variable score_MAX, will store actual biggest score.
						- If compared alignment will have lesser score than score_MAX, variable available will be set to false

						
						ListIterator<Read> iterator2(reads->First());
						while (iterator2.Current() != NULL) {
							Read* rc = iterator2.Current()->Value();
							ListIterator<Alignment> a_iterator_compare(rc->alignments->First());
							while (a_iterator_compare.Current() != NULL) {
								Alignment* a2 = a_iterator_compare.Current()->Value();
								//Start porovnani score
								string x1 = r->sequence;
								string x2 = rc->sequence;
								if (x1 == x2) {
									if (a2->score < this->score_MAX) { //If score is lesser then MAX?score, set available to false
										a2->available = false;
									}
									else if (a2->score == this->score_MAX) {
										//Nothing will happen
									}
									else {
										this->score_MAX = a2->score;
										a->available = false;
									}
								}
								//Konec porovnani score
								a_iterator_compare.Next();
							}

							Read* rc2 = rc->paired_read;
							ListIterator<Alignment> b_iterator_compare(rc2->alignments->First());
							while (b_iterator_compare.Current() != NULL) {
								Alignment* b2 = b_iterator_compare.Current()->Value();
								//Start porovnani score
								string x1 = r->sequence;
								string x2 = rc->sequence;
								if (x1 == x2) {
									if (b2->score < this->score_MAX) { //If score is lesser then MAX?score, set available to false
										b2->available = false;
									}
									else if (b2->score == this->score_MAX) {
										//Nothing will happen
									}
									else {
										this->score_MAX = b2->score;
										a->available = false;
									}
								}
								//Konec porovnani score
								
								b_iterator_compare.Next();
							}
						}
					}
					//End of the control
					*/
					if ((a->score > T) && (a->available == true)) {
						this->ofs << r->name << "\t";
						this->ofs << "FLAG" << "\t";
						this->ofs << a->chromosome << "\t";
						this->ofs << a->pos << "\t";
						this->ofs << "MAPQ" << "\t";
						this->ofs << a->cigar << "\t";
						string x1 = r->name;
						string x2 = r->paired_read->name;
						if (x1 == x2) {
							this->ofs << "=" << "\t";
						}
						else {
							this->ofs << "*" << "\t";
						}
						/*
						this->ofs << b->pos << "\t";
						if (b->pos >= a->pos) {
							temp = b->pos - a->pos;
							temp += b->cigar_length;
						}
						else {
							temp = a->pos - b->pos;
							temp += a->cigar_length;
						}
						this->ofs << temp << "\t";
						*/
						
						if (b_iterator.Current() != NULL) {
							Alignment* b = b_iterator.Current()->Value();
							this->ofs << b->pos << "\t";
							if (b->pos >= a->pos) {
								temp = b->pos - a->pos;
								temp += b->cigar_length;
							}
							else {
								temp = a->pos - b->pos;
								temp += a->cigar_length;
							}
							this->ofs << temp << "\t";
							b_iterator.Next();
						}
						else {
							this->ofs << "N/A" << "\t";
							this->ofs << "N/A" << "\t";
						}
						
						this->ofs << r->sequence << "\t";
						this->ofs << r->quality << "\t";
						this->ofs << "\n";
						temp = 0;
					}
					a_iterator.Next();
					//b_iterator.Next();
					this->score_MAX = 0;
				}

				a_iterator = r->alignments->First();
				b_iterator = r2->alignments->First();
				//while (a_iterator.Current() != NULL && b_iterator.Current() != NULL) {
				while (b_iterator.Current() != NULL) {
					Alignment* b = b_iterator.Current()->Value();
					//Alignment* a = a_iterator.Current()->Value();
					/*
					//Start of the SCORE control to find best match or equal matches.
					if (b->available == true) {
						this->score_MAX = b->score;

						Here is the compare with all alignments, rules below
						- Compare this alignment with all Read & paired read alignments.
						- Variable score_MAX, will store actual biggest score.
						- If compared alignment will have lesser score than score_MAX, variable available will be set to false

						ListIterator<Read> iterator2(reads->First());
						while (iterator2.Current() != NULL) {
							Read* rc = iterator2.Current()->Value();
							ListIterator<Alignment> a_iterator_compare(rc->alignments->First());
							while (a_iterator_compare.Current() != NULL) {
								Alignment* a2 = a_iterator_compare.Current()->Value();
								//Start porovnani score
								string x1 = r->sequence;
								string x2 = rc->sequence;
								if (x1 == x2) {
									if (a2->score < this->score_MAX) { //If score is lesser then MAX?score, set available to false
										a2->available = false;
									}
									else if (a2->score == this->score_MAX) {
										//Nothing will happen
									}
									else {
										this->score_MAX = a2->score;
										b->available = false;
									}
								}
								//Konec porovnani score
								a_iterator_compare.Next();
							}

							Read* rc2 = rc->paired_read;
							ListIterator<Alignment> b_iterator_compare(rc2->alignments->First());
							while (b_iterator_compare.Current() != NULL) {
								Alignment* b2 = b_iterator_compare.Current()->Value();
								//Start porovnani score
								string x1 = r->sequence;
								string x2 = rc->sequence;
								if (x1 == x2) {
									if (b2->score < this->score_MAX) { //If score is lesser then MAX?score, set available to false
										b2->available = false;
									}
									else if (b2->score == this->score_MAX) {
										//Nothing will happen
									}
									else {
										this->score_MAX = b2->score;
										b->available = false;
									}
								}
								//Konec porovnani score

								b_iterator_compare.Next();
							}
						}
					}
					//End of the control
				*/
					if ((b->score > T) && (b->available == true)) {
						this->ofs << r2->name << "\t";
						this->ofs << "FLAG" << "\t";
						this->ofs << b->chromosome << "\t";
						this->ofs << b->pos << "\t";
						this->ofs << "MAPQ" << "\t";
						this->ofs << b->cigar << "\t";
						string x1 = r2->name;
						string x2 = r2->paired_read->name;
						if (x1 == x2) {
							this->ofs << "=" << "\t";
						}
						else {
							this->ofs << "*" << "\t";
						}
						/* second alignment postion (old version)
						this->ofs << a->pos << "\t";
						if (a->pos >= b->pos) {
							temp = a->pos - b->pos;
							temp += a->cigar_length;

						}
						else {
							temp = b->pos - a->pos;
							temp += b->cigar_length;
						}
						this->ofs << -temp << "\t";
						*/

						if (a_iterator.Current() != NULL) {
							Alignment* a = a_iterator.Current()->Value();
							this->ofs << a->pos << "\t";
							if (a->pos >= b->pos) {
								temp = a->pos - b->pos;
								temp += a->cigar_length;

							}
							else {
								temp = b->pos - a->pos;
								temp += b->cigar_length;
							}
							this->ofs << -temp << "\t";
							a_iterator.Next();
						}
						else {
							this->ofs << "N/A" << "\t";
							this->ofs << "N/A" << "\t";
						}
						this->ofs << r2->sequence << "\t";
						this->ofs << r2->quality << "\t";
						this->ofs << "\n";
						temp = 0;
					}
					b_iterator.Next();
					//a_iterator.Next();
					this->score_MAX = 0;
				}

			}
			
			iterator.Next();
		}
		this->ofs.close();
	}
	
};
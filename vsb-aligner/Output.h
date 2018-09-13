#pragma once
#include <fstream>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>

#include "Alignment.h"
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
public:
	Output(char* output_file){
		this->sam_file = output_file;
		int file_check = remove(output_file); //If file exists, return code is 0. We have it here to be sure, that file will be deleted.
	}
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

	void print_program_info(char* ID, char* PN, char* VN, int T, char* r1, char* r2, char* genome) {
		this->ofs.open(this->sam_file, ios::out | ios::app);
		this->ofs << "@PG ID:" << ID << "\tPN:" << PN << "\tVN:" << VN << "\t-T " << T << "\t" << r1 << "\t" << r2 << "\t" << genome << "\n";
		this->ofs.close();
	}

	void print_output_data(List<Read>* reads, int T) {
		this->ofs.open(this->sam_file, ios::out | ios::app);
		ListIterator<Read> iterator(reads->First());
		while (iterator.Current() != NULL) {
			Read* r = iterator.Current()->Value();
			ListIterator<Alignment> a_iterator(r->alignments->First());

			//Version only with paired read - for this work all of the reads is in the pair
			if (r->paired_read != NULL) {
				int temp = 0;
				Read* r2 = r->paired_read;
				ListIterator<Alignment> b_iterator(r2->alignments->First());

				//while (a_iterator.Current() != NULL && b_iterator.Current() != NULL) {
				while (a_iterator.Current() != NULL) {
					Alignment* a = a_iterator.Current()->Value();
					//Alignment* b = b_iterator.Current()->Value();
					if ((a->score) > T) {
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
				}

				a_iterator = r->alignments->First();
				b_iterator = r2->alignments->First();
				//while (a_iterator.Current() != NULL && b_iterator.Current() != NULL) {
				while (b_iterator.Current() != NULL) {
					Alignment* b = b_iterator.Current()->Value();
					//Alignment* a = a_iterator.Current()->Value();
					if ((b->score) > T) {
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
						/*
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
				}

			}
			
			iterator.Next();
		}
		this->ofs.close();
	}
	
};
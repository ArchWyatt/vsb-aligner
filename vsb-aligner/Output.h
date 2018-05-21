#pragma once
#include <fstream>
#include <cstdio>
#include <iostream>

#include "Alignment.h"
#include "Definitions.h"
#include "Genome.h"
#include "GenomeRegion.h"
#include "List.h"

using namespace std;

class Output {
private:

	char* sam_file = NULL;
public:
	Output(char* output_file){
		ofstream ofs(output_file);

		this->sam_file = output_file;
		int file_check = remove(output_file); //If file exists, return code is 0. We have it here to be sure, that file will be deleted.
	}
	void print_head(List<GenomeRegion>* gen){
		ListIterator<GenomeRegion> iterator(gen->First());
		ofstream ofs;
		ofs.open(this->sam_file);
		while (iterator.Current() != NULL) {
			GenomeRegion* temp = iterator.Current()->Value();
			ofs << "@SR SN:" << temp->chromosome_id << " LN:" << temp->bases_number << "\n";
			iterator.Next();
		}
		ofs.close();
	}

	void print_program_info() {
		ofstream ofs;
		ofs.open(this->sam_file);
		// program info
		ofs.close();
	}

	void print_output_data() {

	}
	
};
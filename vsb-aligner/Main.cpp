#include <iostream>
#include <cstring>
#include <chrono>

#include "Aligner.h"
#include "Definitions.h"
#include "Genome.h"
#include "Read.h"
#include "Needleman_Wunch.h"
#include "Needleman_Wunch_alternative.h"
#include "Smith_Waterman.h"
#include "Output.h"
#include "MAPQ.h"

using namespace std;

void PrintUsage()
{
	cerr << "Usage:" << endl << endl;
	cerr << "Aligner -r1 [forward reads] -r2 [reverse reads] -g [genome_file] -s [SAM_file]" << endl << endl;
	cerr << "Options:" << endl;
	cerr << "\t -r1 -forward reads *_R1_* file; Mandatory." << endl;
	cerr << "\t -r2 -forward reads *_R2_* file; Mandatory." << endl;
	cerr << "\t -g -input genome; Default: None" << endl;	
	cerr << "\t -s -output name of SAM file; Default: None" << endl;
}

void ParseArguments(int argc, char* argv[], ProgInfo* prog_info)
{
	u_int arg_id = 1;
	while (arg_id < argc) {
		if (strcmp(argv[arg_id], "-r1") == 0) {
			prog_info->fq_F = argv[arg_id + 1];			
		}
		else if (strcmp(argv[arg_id], "-r2") == 0) {
			prog_info->fq_R = argv[arg_id + 1];			
		}
		else if (strcmp(argv[arg_id], "-g") == 0) {
			prog_info->genome_path = argv[arg_id + 1];			
		}
		else if (strcmp(argv[arg_id], "-s") == 0) {
			prog_info->sam_file = argv[arg_id + 1];
		}
		else {
			cerr << "Unknown option: " << argv[arg_id] << endl;
			PrintUsage();
			exit(EXIT_FAILURE);
		}
		arg_id += 2;
	}
}

void PrintArgs(ProgInfo* prog_info)
{
	cout << "Aligner: vsb-aligner" << endl;
	cout << "Reads R1: " << prog_info->fq_F << endl;
	cout << "Reads R2: " << prog_info->fq_R << endl;
	cout << "Genome: " << prog_info->genome_path << endl;
	cout << "SAM file: " << prog_info->sam_file << endl << endl;
}

void Verify(ProgInfo* proginfo)
{
	if (proginfo->genome_path == NULL) {
		cerr << "Genome path is missing." << endl;
		PrintUsage();
		exit(EXIT_FAILURE);
	}
	else if (proginfo->fq_F == NULL) {
		cerr << "Reads R1 are missing." << endl;
		PrintUsage();
		exit(EXIT_FAILURE);
	}
	else if (proginfo->fq_R == NULL) {
		cerr << "Reads R2 are missing." << endl;
		PrintUsage();
		exit(EXIT_FAILURE);
	}
	else if (proginfo->sam_file == NULL) {
		cerr << "SAM file name missing." << endl;
		PrintUsage();
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[])
{
	/* class holding program specific informations */
	ProgInfo prog_info;

	/* parse arguments */
	ParseArguments(argc, argv, &prog_info);
		
	/* verify arguments */
	Verify(&prog_info);

	/* print arguments */
	PrintArgs(&prog_info);
	
	/* verify existence of genome and fai  file */
	if (!Genome::GenomeExists(prog_info.genome_path)) {
		cerr << "Genome path: " << prog_info.genome_path << " doesn't exist." << endl;
		exit(EXIT_FAILURE);
	}

	if (!Genome::GenomeIndexExists(prog_info.genome_path)) {
		cerr << "Index for genome: " << prog_info.genome_path << " doesn't exist in the same folder." << endl;
		exit(EXIT_FAILURE);
	}
	
	Genome genome(prog_info.genome_path);
	
	//zkontrolovat existenci sa-indexu
	if (!Genome::SAIndexExists(prog_info.genome_path)) {
		cout << "Suffix array not found." << endl << "Beginning suffix array build... :" << endl;
		genome.PrepareIndexes();
	}
	
	//check indexes
	genome.CheckSAIndexes();
	
	Aligner aligner(&prog_info, &genome);
	
	//will do the pairing of reads
	aligner.PairReads();
	
	//will assign each read (0,1 or multiple) alignments
	aligner.AlignReads();
	
	//nacist do pameti ready
	List<Read>* reads = aligner.Reads();
	
	/*
	Sekce: Martin Kubala
	*/
	//Start time measuring
	auto start = chrono::high_resolution_clock::now();
	cout << "Martin Kubala - Part" << endl;

	//SAM file header
	Output *output = new Output(prog_info.sam_file);
	List<GenomeRegion>* genom_out = genome.Chromosomes();
	output->print_head(genom_out);
	output->print_program_info(prog_info.options->ID, prog_info.options->PN, prog_info.options->VN, prog_info.options->T, prog_info.fq_F, prog_info.fq_R, prog_info.genome_path);

	//Start of the computing part
	cout << "Computing part" << endl;
	
	//Read list iteration
	ListIterator<Read> iterator(reads->First());
	while (iterator.Current() != NULL){
		Read* r = iterator.Current()->Value();

		//Read alignment list iteration
		ListIterator<Alignment> a_iterator(r->alignments->First());
		while (a_iterator.Current() != NULL) {
			Alignment* a = a_iterator.Current()->Value();
			char *reference_genome = genome.BaseIntervalDisc(a->chromosome, ((a->pos) - (prog_info.options->range_prefix)), ((a->pos + r->seq_len - 1) + prog_info.options->range_suffix));

			Smith_Waterman *test = new Smith_Waterman(r->sequence, reference_genome, prog_info.options->gap_score, prog_info.options->match_score, prog_info.options->mismatch_score);
			//Needleman_Wunch *test = new Needleman_Wunch(r->sequence, reference_genome, prog_info.options->gap_score, prog_info.options->match_score, prog_info.options->mismatch_score);
			//Needleman_Wunch_alternative *test = new Needleman_Wunch_alternative(r->sequence, reference_genome, prog_info.options->gap_score, prog_info.options->match_score, prog_info.options->mismatch_score);
			a->pos = a->pos + (test->get_first_pos() - 1);
			a->cigar = test->get_cigar();
			a->cigar_length = test->get_cigar_length();
			a->score = test->get_matrix_max_score();

			//MAP Quality part
			MAPQ *temp_MAPQ1 = new MAPQ(test->get_mismatch(), r->quality);
			a->MAPQ = temp_MAPQ1->get_MAPQ();

			//FLAG part
			a->FLAG += 1;	// + 0x1 - template having multiple segments in sequencing - This part was done in preprocessed data
			a->FLAG += 2;	// + 0x2 - each segment properly aligned according to the aligner - This part was done in preprocessed data
			a->FLAG += 32;	// + 0x20 - SEQ of the next segment in the template being reverse complemented - This part was done in preprocessed data, all reads are paired.
			a->FLAG += 64;	// + 0x40 - the ﬁrst segment in the template

			//Clearing allocated memory
			delete temp_MAPQ1;
			delete test;
			delete[] reference_genome;

			a_iterator.Next();
		}

		//Paired read
		Read* r2 = r->paired_read;
		
		//Paired Read alignment list iteration
		ListIterator<Alignment> b_iterator(r2->alignments->First());
		while (b_iterator.Current() != NULL) {
			Alignment* b = b_iterator.Current()->Value();
			char *reference_genome2 = genome.BaseIntervalDisc(b->chromosome, ((b->pos) - (prog_info.options->range_prefix)), ((b->pos + r2->seq_len - 1) + prog_info.options->range_suffix));
			Smith_Waterman *test2 = new Smith_Waterman(r2->sequence, reference_genome2, prog_info.options->gap_score, prog_info.options->match_score, prog_info.options->mismatch_score);
			//Needleman_Wunch *test2 = new Needleman_Wunch(r2->sequence, reference_genome2, prog_info.options->gap_score, prog_info.options->match_score, prog_info.options->mismatch_score);
			//Needleman_Wunch_alternative *test2 = new Needleman_Wunch_alternative(r2->sequence, reference_genome2, prog_info.options->gap_score, prog_info.options->match_score, prog_info.options->mismatch_score);
			b->pos = b->pos + (test2->get_first_pos() - 1);
			b->cigar = test2->get_cigar();
			b->cigar_length = test2->get_cigar_length();
			b->score = test2->get_matrix_max_score();

			//MAP Quality part
			MAPQ *temp_MAPQ2 = new MAPQ(test2->get_mismatch(), r2->quality);
			b->MAPQ = temp_MAPQ2->get_MAPQ();

			//FLAG part
			b->FLAG += 1;	// + 0x1 - template having multiple segments in sequencing - This part was done in preprocessed data
			b->FLAG += 2;	// + 0x2 - each segment properly aligned according to the aligner - This part was done in preprocessed data
			b->FLAG += 16;	// + 0x10 - SEQ being reverse complemented - This part was done in preprocessed data, all reads are paired.
			b->FLAG += 128;	// + 0x80 - the last segment in the template 

			//Clearing allocated memory
			delete temp_MAPQ2;
			delete test2;
			delete[] reference_genome2;

			b_iterator.Next();
		}
		iterator.Next();
	}
	cout << "SAM output part" << endl;
	output->output_prepare(reads);
	output->output_top_score_filtering(prog_info.options->T);
	cout << "Printing into file" << endl;
	output->print_output_data();

	auto elapsed = chrono::high_resolution_clock::now() - start;
	long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
	long long milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
	long long seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
	long long minutes = std::chrono::duration_cast<std::chrono::minutes>(elapsed).count();
	cout << "Runtime in microseconds: " << microseconds << endl;
	cout << "Runtime in miliseconds: " << milliseconds << endl;
	cout << "Runtime in seconds: " << seconds << endl;
	cout << "Runtime in minutes: " << minutes << endl;
	system("pause");

	return EXIT_SUCCESS;
}
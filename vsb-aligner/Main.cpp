#include <iostream>
#include <cstring>
#include <chrono>

#include "Aligner.h"
#include "Definitions.h"
#include "Genome.h"
#include "Read.h"
#include "Needleman_Wunch.h"
#include "Needleman_Wunch_Old.h"
#include "Smith_Waterman.h"

using namespace std;

void PrintUsage()
{
	cerr << "Usage:" << endl << endl;
	cerr << "Aligner -r1 [forward reads] -r2 [reverse reads] -g [genome_file]" << endl << endl;
	cerr << "Options:" << endl;
	cerr << "\t -r1 -forward reads *_R1_* file; Mandatory." << endl;
	cerr << "\t -r2 -forward reads *_R2_* file; Mandatory." << endl;
	cerr << "\t -g -input genome; Default: None" << endl;	
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
	cout << "Genome: " << prog_info->genome_path << endl << endl;	
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
	
	//zkotnrolovat existenci sa-indexu
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
	cout << "Zde zacina ma cast" << endl;

	//system("pause");
	/*
	Sekce: Martin Kubala
	*/

	auto start = chrono::high_resolution_clock::now();
	u_int rozsah = 20; //Rozsah genom zvetsen o n znaku pred a n znaku po genomu pro vetsi presnost urceni pozice genomu.

	int gap_score = -2;
	int match_score = 5;
	int mismatch_score = -1;

	//Needleman_Wunch *test = new Needleman_Wunch(a, b, gap_score, match_score, mismatch_score);
	//Needleman_Wunch_Old *test = new Needleman_Wunch_Old(a, b, gap_score, match_score, mismatch_score);
	
	ListIterator<Read> iterator(reads->First());
	while (iterator.Current() != NULL){
		Read* r = iterator.Current()->Value();

		ListIterator<Alignment> a_iterator(r->alignments->First());

		while(a_iterator.Current() != NULL){
			Alignment* a = a_iterator.Current()->Value();
			//cigar ulozit do aligmentu

			//cout << a->chromosome << "\t" << a->pos << endl;
			char *reference_genome = genome.BaseIntervalDisc(a->chromosome, (a->pos - rozsah), ((a->pos + r->seq_len - 1) + rozsah));   //Zde bych přidal ten interval +/- 20? reference genome OK kdyz dam interval +/- 20
									  
			Smith_Waterman *test = new Smith_Waterman(r->sequence, reference_genome, gap_score, match_score, mismatch_score);

			system("pause");

			Read* r2 = r->paired_read;	
			Smith_Waterman *test2 = new Smith_Waterman(r2->sequence, reference_genome, gap_score, match_score, mismatch_score);

			system("pause");
			a_iterator.Next();
		}

		iterator.Next();
	}

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
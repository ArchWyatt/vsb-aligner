#include <iostream>
#include <cstring>

#include "Definitions.h"
#include "Genome.h"
#include "SuffixArray.h"

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
	cout << "Aligner: " << endl;
	cout << "Reads R1: " << prog_info->fq_F;
	cout << "Reads R2: " << prog_info->fq_R;
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
	//ParseArguments(argc, argv, &prog_info);

	/* verify arguments */
	//Verify(&prog_info);

	/* verify existence of genome and fai  file */
	//if (!Genome::GenomeExists(prog_info.genome_path)) {
	//	cerr << "Genome path: " << prog_info.genome_path << " doesn't exist." << endl;
	//	exit(EXIT_FAILURE);
	//}
	//if (!Genome::GenomeIndexExists(prog_info.genome_path)) {
	//	cerr << "Index for genome: " << prog_info.genome_path << " doesn't exist in the same folder." << endl;
	//	exit(EXIT_FAILURE);
	//}

	//Genome genome(prog_info.genome_path);
	
	char* test = "TTTT";
	SuffixArray sa(test, strlen(test)+1);	
	sa.Verify();

	return EXIT_SUCCESS;
}
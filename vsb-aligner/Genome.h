#pragma once
#ifndef __GENOME_H__
#define __GENOME_H__

/*
	Represents operations with genome.

	Functions:
		- provides acces to proper positions of chromosomes - will be searched sequentially
		- prefetch chromosome into memory

	TODO:
		- fasta fai file should not be const
		- delete lists
		- reverse complement - alphabet should be accessed in constant way, not by if-ing
*/

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>

#include "Definitions.h"
#include "GenomeRegion.h"
#include "Hash.h"
#include "List.h"
#include "Utils.h"

using namespace std;


class Genome
{
private:
	const char* GENOME_INDEX_SUFFIX = ".fai";

	char* genome_path;
	char* index_path;

	/*
		List of chromosomes - positions and length, corresponds to fasta.fai
	*/
	List<GenomeRegion>* chromosome_list;
	Hash<GenomeRegion>* active_chromosomes_hash;

	/* 
		Loads bases for particular chromosome into memory.
	*/
	void FetchChromosome(GenomeRegion* region);

public:
	Genome(char* path_to_genome);
	~Genome();

	/*
		Initialize Genome List
	*/
	void InitializeGenomeList();

	/* 
		Prefetches chromosomes bases into memory.
	*/
	void PrefetchChromosomes(List<char>* active_chromosomes);

	/*
		Prefetech one chromosome.
	*/
	void PrefetchChromosome(char* chromosome);

	/* 
		Returns base in chrX:pos system.
	*/
	char Base(char* chr, u_int pos);

	/* 
		Returns interval of bases between(including) two positions.
	*/
	char* BaseInterval(char* chr, u_int pos_start, u_int pos_end);

	/*
		Returns interval with no prefetch of chromosome, uses direct load from disk.
	*/
	char* BaseIntervalDisc(char* chr, u_int pos_start, u_int pos_end);

	/*
		Returns the list of chromosomes.
	*/
	List<GenomeRegion>* Chromosomes();

	/*
		Returns reverse complement of the sequence. Operates on itself
	*/
	static char ReverseComplement(char base);

	/*
		Returns reverse complement of the sequence. Operates and returns itself.
	*/
	static char* ReverseComplement(char* sequence);

	/*
		Returns new reverse complement of the sequence.
	*/
	static char* NewReverseComplement(char* sequence);
		
	/*
		Verifies genome existence.
	*/
	static bool GenomeExists(char* genome_path);

	/*
		Verifies fasta.fai index existence.
	*/
	static bool GenomeIndexExists(char* genome_path);
};
#endif
#include "Genome.h"

Genome::Genome(char* path_to_genome)
{
	genome_path = path_to_genome;
	
	index_path = Utils::StrAppend(genome_path, GENOME_INDEX_SUFFIX);

	chromosome_list = new List<GenomeRegion>();
	
	InitializeGenomeList();
}

Genome::~Genome()
{
	ListIterator<GenomeRegion> list_iterator(chromosome_list->First());
	while (list_iterator.Current() != NULL) {
		delete list_iterator.Current()->Value();
		list_iterator.Next();
	}

	delete chromosome_list;
	delete active_chromosomes_hash;
}

void Genome::InitializeGenomeList()
{
	ifstream chrom_file(index_path, ios::binary);
	
	char line[1024];

	while (!chrom_file.eof()) {
		chrom_file.getline(line,1024,'\n');
		chromosome_list->Append(Utils::LineToRegion(line));
	}

	chrom_file.close();
}

void Genome::FetchChromosome(GenomeRegion* region)
{
	//open the genome file
	ifstream genome(genome_path, ios::binary);
	if (genome.fail()) {
		cerr << "ERROR: cannot read genome file." << endl;
		exit(EXIT_FAILURE);
	}

	//prepare the size of the chrom buffer
	//further bases will c-string zero terminated
	region->bases = new char[region->bases_number + 1];
	//make the zero-th element be max byte value	
	region->bases[region->bases_number] = 0;

	//set the file ptr to the location given by genome region object
	genome.seekg(region->bases_start);
	
	u_int lines_no = region->bases_number / region->bases_per_line;	
	u_int rest = region->bases_number % region->bases_per_line;

	u_int region_ptr = 0;
	u_int line_ptr = 0;
	
	//read line by line	and fetch into memory	
	char* line = new char[region->bytes_per_line];
	while (line_ptr < lines_no) {		
		genome.read(line, region->bytes_per_line);		
		Utils::StrUppercase(line, region->bytes_per_line);
		//place line into region buffer		
		memcpy(region->bases + line_ptr*region->bases_per_line, line, region->bases_per_line * sizeof(char));
		//skip the end of line		
		line_ptr++;
	}

	if (rest > 0) {
		//read the rest of bases		
		genome.read(line, rest);		
		Utils::StrUppercase(line, region->bytes_per_line);		
		memcpy(region->bases + line_ptr*region->bases_per_line, line, rest * sizeof(char));
	}

	genome.close();
}

void Genome::ReleaseChromosome(GenomeRegion* region)
{
	/*
		If suffix array already built, then we have to delete bases manualy.
	*/
	//if(suffix_built)
	delete[] region->bases;

	region->bases = NULL;
}

void Genome::PrefetchChromosome(char* chromosome)
{
	ListIterator<GenomeRegion> region_iterator(chromosome_list->First());

	while (region_iterator.Current() != NULL) {
		GenomeRegion* region = region_iterator.Current()->Value();

		if (strcmp(region->chromosome_id, chromosome) == 0) {
			FetchChromosome(region);
			active_chromosomes_hash->Add(chromosome, region);
			break;
		}

		region_iterator.Next();
	}
}

void Genome::PrefetchChromosomes(List<char>* active_chromosomes)
{
	//find particular chromosome in the chromosome list
	ListIterator<GenomeRegion> regions_iterator(chromosome_list->First());
	ListIterator<char> active_iterator(active_chromosomes->First());

	while (active_iterator.Current() != NULL) {
		while (regions_iterator.Current() != NULL) {
			char* active_region_name = active_iterator.Current()->Value();
			char* region_name = regions_iterator.Current()->Value()->chromosome_id;

			if (strcmp(active_region_name, region_name) == 0) {
				//match => FetchChromosome from genome file
				FetchChromosome(regions_iterator.Current()->Value());
				//place region into the active chromosomes hash
				active_chromosomes_hash->Add(region_name, regions_iterator.Current()->Value());				
				break;
			}

			regions_iterator.Next();
		}

		active_iterator.Next();
	}
}

char Genome::Base(char* chr, u_int pos)
{
	GenomeRegion* region = active_chromosomes_hash->Get(chr);

	if (region == NULL) {
		PrefetchChromosome(chr);
		region = active_chromosomes_hash->Get(chr);
	}

	if (pos >= region->bases_number) {
		cerr << "ERROR: demanding position number large than the size of the chromosome." << endl;
		exit(EXIT_FAILURE);
	}

	return region->bases[pos];
}

char* Genome::BaseInterval(char* chr, u_int pos_start, u_int pos_end)
{
	GenomeRegion* region = active_chromosomes_hash->Get(chr);

	if (region == NULL){
		PrefetchChromosome(chr);
		region = active_chromosomes_hash->Get(chr);
	}

	if (pos_end < pos_start) {
		cerr << "ERROR: position of the beginning of the interval is greater than the position of the end." << endl;
		exit(EXIT_FAILURE);
	}

	if (pos_end >= region->bases_number) {
		cerr << "ERROR: demanding interval end number is larger than the size of the chromosome." << endl;
		exit(EXIT_FAILURE);
	}

	//+1 beacuse the last position is included
	u_int i_len = pos_end - pos_start + 1;
	//allocate space
	char* interval = new char[i_len + 1];
	//make it a C-string
	interval[i_len] = 0;
	memcpy(interval, region->bases + pos_start, sizeof(char)*i_len);

	return interval;
}

char* Genome::BaseIntervalDisc(char* chromosome, u_int pos_start, u_int pos_end)
{
	GenomeRegion* chrom = NULL;

	// find genomic region
	ListIterator<GenomeRegion> region_iterator(chromosome_list->First());

	while (region_iterator.Current() != NULL) {
		GenomeRegion* region = region_iterator.Current()->Value();

		if (strcmp(region->chromosome_id, chromosome) == 0) {
			chrom = region;
			break;
		}

		region_iterator.Next();
	}

	// compute beginning of the interval in the chromosome
	u_int start_line = pos_start / chrom->bases_per_line;
	u_int pos_in_line = pos_start % chrom->bases_per_line;

	u_int pos_in_file = chrom->bases_start + start_line * chrom->bytes_per_line + pos_in_line - 1;

	u_int buf_len = pos_end - pos_start + 1;
	char* outbuf = new char[buf_len + 1];
	outbuf[buf_len] = 0;

	ifstream input(genome_path, ios::binary);
	input.seekg(pos_in_file);

	u_int buf_ptr = 0;
	while (buf_ptr < buf_len) {
		//compute the number of bases to read to the buffer
		if (pos_in_line != 0 && chrom->bases_per_line - pos_in_line >= buf_len) {
			input.read(outbuf, sizeof(char) * (buf_len));
			buf_ptr += chrom->bases_per_line - pos_in_line + 1;
			pos_in_file = input.tellg();
			pos_in_file++;
			input.seekg(pos_in_file);
			pos_in_line = 0;
		}
		else if (pos_in_line != 0) {
			input.read(outbuf, sizeof(char) * (chrom->bases_per_line - pos_in_line + 1));
			buf_ptr += chrom->bases_per_line - pos_in_line + 1;
			pos_in_file = input.tellg();
			pos_in_file++;
			input.seekg(pos_in_file);
			pos_in_line = 0;
		}
		else if (buf_len - buf_ptr < chrom->bases_per_line) {
			input.read(outbuf + buf_ptr, sizeof(char) * (buf_len - buf_ptr));
			buf_ptr += buf_len - buf_ptr;
		}
		else {
			input.read(outbuf + buf_ptr, sizeof(char) * chrom->bases_per_line);
			buf_ptr += chrom->bases_per_line;
			pos_in_file = input.tellg();
			pos_in_file++;
			input.seekg(pos_in_file);
			pos_in_line = 0;
		}
	}

	Utils::StrUppercase(outbuf, buf_len);

	return outbuf;
}

List<GenomeRegion>* Genome::Chromosomes()
{
	return chromosome_list;
}

void Genome::PrepareIndexes()
{
	char* sa_index_name = Utils::StrAppend(this->genome_path, ".sa");
	char* bwt_index_name = Utils::StrAppend(this->genome_path, ".bwt");

	ListIterator<GenomeRegion> iterator(chromosome_list->First());

	while (iterator.Current() != NULL) {
		GenomeRegion* region = iterator.Current()->Value();
		cout << region->chromosome_id << endl;

		FetchChromosome(region);

		SuffixArray sa(region->bases, region->bases_number + 1);
		sa.Verify();
		sa.Serialize(sa_index_name);
		//sa.SerializeBWT(bwt_index_name);

		ReleaseChromosome(region);
		iterator.Next();
	}
}

void Genome::CheckSAIndexes()
{
	char* sa_index_name = Utils::StrAppend(this->genome_path, ".sa");
	ifstream sa_file(sa_index_name, ios::binary);

	ListIterator<GenomeRegion> iterator(chromosome_list->First());
	long long start_pos = 0;
	while (iterator.Current() != NULL) {
		GenomeRegion* region = iterator.Current()->Value();

		u_int* test = new u_int[1];

		sa_file.seekg(start_pos);

		sa_file.read((char*)test, 4);

		if(test[0] != region->bases_number){
			cout << region->chromosome_id << endl;
			cout << "Err: " << test[0] << " \t " << region->bases_number << endl;
		}

		start_pos += 4 * (region->bases_number + 1);
		delete[] test;

		iterator.Next();
	}
}

char Genome::ReverseComplement(char base)
{
	if (base == 'A')
		return 'T';
	else if (base == 'C')
		return 'G';
	else if (base == 'G')
		return 'C';
	else if (base == 'T')
		return 'A';
	else {
		cerr << "ERROR: Symbol " << base << " is not base symbol." << endl;
		exit(EXIT_FAILURE);
	}
}

char* Genome::ReverseComplement(char* sequence)
{
	u_int s_len = strlen(sequence);
	//reverse symbols
	for (u_int i = 0; i < s_len / 2; i++) {
		char aux = sequence[s_len - i - 1];
		sequence[s_len - i - 1] = sequence[i];
		sequence[i] = aux;
	}
	//each symbol complemented
	for (u_int i = 0; i < s_len; i++)
		sequence[i] = ReverseComplement(sequence[i]);

	return sequence;
}

char* Genome::NewReverseComplement(char* sequence)
{
	u_int s_len = strlen(sequence);
	char* reversed_sequence = new char[s_len + 1];
	//c-string
	reversed_sequence[s_len] = 0;
	//reverse symbols
	for (u_int i = 0; i < s_len; i++) {
		reversed_sequence[s_len - 1 - i] = sequence[i];
	}

	//each symbol complemented
	for (u_int i = 0; i < s_len; i++)
		reversed_sequence[i] = ReverseComplement(reversed_sequence[i]);

	return reversed_sequence;
}

bool Genome::GenomeExists(char* genome_path)
{
	ifstream input(genome_path);
	if(input.fail()){	
		input.close();
		return false;
	}
	else {
		input.close();
		return true;
	}
}

bool Genome::GenomeIndexExists(char* genome_path)
{
	// convert genome path to fasta fai path
	const char* index_path = Utils::StrAppend(genome_path, ".fai");

	ifstream input(index_path);
	if (input.fail()) {
		input.close();
		return false;
	}
	else {
		input.close();
		return true;
	}
}

bool Genome::SAIndexExists(char* genome_path)
{
	// convert genome path to sa path
	const char* index_path = Utils::StrAppend(genome_path, ".sa");

	ifstream input(index_path);
	if (input.fail()) {
		input.close();
		return false;
	}
	else {
		input.close();
		return true;
	}
}
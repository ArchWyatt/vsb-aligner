#include "Aligner.h"

Aligner::Aligner()
{
}

Aligner::Aligner(ProgInfo* prog_info, Genome* genome) : m_prog_info(prog_info), m_genome(genome)
{
	forward_reads = new List<Read>();
	reverse_reads = new List<Read>();

	//InitializeReads(prog_info->fq_F, FORWARD_READ);
	//InitializeReads(prog_info->fq_R, REVERSE_READ);
}

Aligner::~Aligner()
{
	/*
	ListIterator<Read> list_iterator_f(forward_reads->First());

	while (list_iterator_f.Current() != NULL) {
		delete list_iterator_f.Current()->Value();
		list_iterator_f.Next();
	}

	delete forward_reads;

	ListIterator<Read> list_iterator_r(reverse_reads->First());

	while (list_iterator_r.Current() != NULL) {
		delete list_iterator_r.Current()->Value();
		list_iterator_r.Next();
	}

	delete reverse_reads;
	*/
}

void Aligner::InitializeReads(char* reads_path, u_char read_type)
{
	ifstream reads(reads_path, ios::binary);

	if (reads.fail()) {
		cerr << "ERROR: cannot open reads fastq file." << endl;
		exit(EXIT_FAILURE);
	}

	char line[1024];
	u_int line_len;
	u_int i = 0;
	//first read into line
	reads.getline(line, 1024);
	while (!reads.eof()) {
		//read four lines related to each read
		Read* read = new Read();

		//read NGS marking string
		line_len = strlen(line);
		if (line[line_len - 1] == '\r')
			line_len--;

		char* ngs_identifier = new char[line_len + 1];
		ngs_identifier[line_len] = 0;
		memcpy(ngs_identifier, line, sizeof(char)*line_len);
		read->descriptor = ngs_identifier;

		//read bases
		reads.getline(line, 1024);
		line_len = strlen(line);
		if (line[line_len - 1] == '\r')
			line_len--;
		char* read_bases = new char[line_len + 1];
		read_bases[line_len] = 0;
		memcpy(read_bases, line, sizeof(char)*line_len);
		read->sequence = read_bases;
		read->seq_len = strlen(read->sequence);

		//read comment
		reads.getline(line, 1024);
		line_len = strlen(line);
		if (line[line_len - 1] == '\r')
			line_len--;
		char* read_comment = new char[line_len + 1];
		read_comment[line_len] = 0;
		memcpy(read_comment, line, sizeof(char)*line_len);
		read->comment = read_comment;

		//read qualitites
		reads.getline(line, 1024);
		line_len = strlen(line);
		if (line[line_len - 1] == '\r')
			line_len--;
		char* read_qualities = new char[line_len + 1];
		read_qualities[line_len] = 0;
		memcpy(read_qualities, line, sizeof(char)*line_len);
		read->quality = read_qualities;
		
		read->ParseIDs();

		if (read_type == FORWARD_READ) {
			read->read_type = FORWARD_READ;
			forward_reads->Append(read);
		}
		else if (read_type == REVERSE_READ) {
			read_type = REVERSE_READ;		
			reverse_reads->Append(read);
		}
		i++;
		//the last read of the identifier => if not loaded => eof
		reads.getline(line, 1024);
	}

	reads.close();
}

void Aligner::PairReads()
{
	ListIterator<Read> forward_iterator(forward_reads->First());
	ListIterator<Read> reverse_iterator(reverse_reads->First());

	while (forward_iterator.Current() != NULL && reverse_iterator.Current() != NULL) {
		Read* f_read = forward_iterator.Current()->Value();
		Read* r_read = reverse_iterator.Current()->Value();
		
		if (f_read->read_idx1 == r_read->read_idx1 && f_read->read_idx2 == r_read->read_idx2) {
			f_read->paired_read = r_read;
			r_read->paired_read = f_read;				
		}		
		else {
			cerr << "Error: Aligner[PairReads] - reads not paired." << endl;
			exit(EXIT_FAILURE);
		}

		reverse_iterator.Next();
		forward_iterator.Next();
	}
}

void Aligner::AlignReads()
{
	/*
		0) pass each chromosome
		1) pass each read
		2) try to localize full length sequence using sa and bwt
			found 1 or more => then append alignment(s) => next read
			not found => step 3.
		3) pigeonhole => candidate alignments => given the initial setup min_number of matches 3/|splits| and split size 20
	*/

	List<GenomeRegion>* chromosomes = m_genome->Chromosomes();
	ListIterator<GenomeRegion> chromosome_iterator(chromosomes->First());

	char* sa_path = Utils::StrAppend(m_prog_info->genome_path, ".sa");
	
	long long chrom_start = 0;
	while (chromosome_iterator.Current() != NULL) {
		GenomeRegion* chromosome = chromosome_iterator.Current()->Value();
		
		cout << "Aligning in:" << chromosome->chromosome_id << endl;
					
		m_genome->FetchChromosome(chromosome);
		
		//load suffix array and prepare bwt		
		SuffixArray sa(sa_path, chromosome->bases, chrom_start, chromosome->bases_number + 1);

		ListIterator<Read> reads_iterator(forward_reads->First());
		while (false && reads_iterator.Current() != NULL) {
			Read* r = reads_iterator.Current()->Value();
			u_int occ = 0;
			u_int* positions = sa.Localize(r->sequence, strlen(r->sequence), occ);

			if (positions != NULL) {
				for (u_int i = 0; i < occ; i++)
					r->alignments->Append(new Alignment(chromosome->chromosome_id, positions[i]));
			}
			else {
				//pigeonhole principle
				//find the minimum size split larger than min_split
				u_char groups = 1;
				while (r->seq_len / (groups + 1) > m_prog_info->aligner_min_split_size)
					groups++;

				u_char split_size = r->seq_len / groups;

				//we know the number of groups
				u_int** group_positions = new u_int*[groups];
				u_int* group_len = new u_int[groups];

				u_int seq_start = 0;
				u_int tot_rec = 0;
				for (u_int i = 0; i < groups; i++) {
					u_int occ = 0;
					u_int* positions = NULL;
					if (i == groups - 1)
						positions = sa.Localize(r->sequence + seq_start, r->seq_len - seq_start, occ);
					else 
						positions = sa.Localize(r->sequence + seq_start, split_size, occ);

					group_positions[i] = positions;
					group_len[i] = occ;
					tot_rec += occ;
					seq_start += split_size;
				}

				//merge all positions into one array
				if (tot_rec > 0){
					u_int* all_positions = new u_int[tot_rec];
					u_int cur_pos = 0;
					for (u_int i = 0; i < groups; i++) {
						for (u_int j = 0; j < group_len[i]; j++)
							all_positions[cur_pos++] = group_positions[i][j];
					}

					qsort(all_positions, tot_rec, sizeof(u_int), Compare);

					u_int search_id = 1;
					u_int current = all_positions[0];
					u_char current_length = 1;
					while (search_id != tot_rec) {
						if (all_positions[search_id] == current)
							current_length++;						
						else {
							//check the current length and potentially append new candidate alignment
							if (current_length > groups / 2)
								r->alignments->Append(new Alignment(chromosome->chromosome_id, current));							
							//reset
							current = all_positions[search_id];
							current_length = 1;
						}
						search_id++;
					}
					delete[] all_positions;
				}

				for (u_int i = 0; i < groups; i++)
					delete[] group_positions[i];

				delete[] group_positions;
				delete[] group_len;
			}
			reads_iterator.Next();
		}

		chrom_start += chromosome->bases_number + 1;
		chromosome_iterator.Next();
	}
}

List<Read>* Aligner::Reads()
{
	return forward_reads;
}

int Aligner::Compare(const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}
#include "SuffixArray.h"

SuffixArray::SuffixArray()
{
	m_rank = NULL;
	m_offsets = NULL;
	m_bwt = NULL;
}

SuffixArray::SuffixArray(char* sequence, u_int sequence_length) : SuffixArray()
{
	m_sequence = sequence;
	m_sequence_length = sequence_length;

	s_array = new u_int[sequence_length];
	for (u_int i = 0; i < sequence_length; i++)
		s_array[i] = i;

	baseMap = new char[256];
	baseMap[0] = 0;
	baseMap['A'] = 1;
	baseMap['C'] = 2;
	baseMap['G'] = 3;
	baseMap['N'] = 4;
	baseMap['T'] = 5;

	SortArray(0, m_sequence_length - 1);
}

SuffixArray::SuffixArray(char* index_path, char* bases, long long pos_start, u_int length) : SuffixArray()
{
	m_sequence = bases;
	m_sequence_length = length;
	s_array = new u_int[m_sequence_length];	
	
	//load suffix array from file
	ifstream input(index_path, ios::binary);
	input.seekg(pos_start*sizeof(u_int));
	input.read((char*)s_array, sizeof(u_int)*length);
	input.close();
	
	PrepareBWT();
}

SuffixArray::~SuffixArray()
{
	delete[] s_array;
	//delete[] m_sequence;
	delete[] m_bwt;
	delete[] m_offsets;
	delete[] m_rank;	
}

void SuffixArray::PrepareBWT()
{
	m_bwt = new char[m_sequence_length];
	m_rank = new u_int[m_sequence_length];
	m_offsets = new u_int[256];

	for (u_int i = 0; i < m_sequence_length; i++) {
		if (s_array[i] == 0)
			m_bwt[i] = 0;
		else
			m_bwt[i] = m_sequence[s_array[i] - 1];
	}

	u_int f[256];
	for (u_int i = 0; i < 256; i++){
		m_offsets[i] = 0;
		f[i] = 0;
	}
	
	for (u_int i = 0; i < m_sequence_length; i++){
		m_rank[i] = f[m_bwt[i]];
		f[m_bwt[i]]++;
	}
		
	m_offsets[0] = 0;
	for (u_int i = 1; i < 256; i++)
		m_offsets[i] = m_offsets[i - 1] + f[i-1];
	
	//delete[] m_sequence;
	//m_sequence = NULL;
}

void SuffixArray::RadixSortRecursive(u_int start, u_int end, u_int level)
{
	//prepare frequency array
	u_int f[SA_ALFA_LEN];
	for (u_int i = 0; i < SA_ALFA_LEN; i++)
		f[i] = 0;

	//count frequencies in interval
	for (u_int i = start; i <= end; i++)
		f[baseMap[m_sequence[s_array[i]+level]]]++;

	//offset array - pointers to the beginning of the sa interval for particular symbol
	u_int off[SA_ALFA_LEN];
	u_int off_end[SA_ALFA_LEN];
	off[0] = start;
	off_end[0] = start + f[0];

	u_int tot = start + f[0];
	for (u_int i = 1; i < SA_ALFA_LEN; i++) {
		off[i] = tot;
		off_end[i] = off[i] + f[i];
		tot += f[i];
	}

	/*
		Postup:
		1 - prochazime po jednotlivych offsetech
		2 - umistujeme prvky na spravnou pozici, dokud nenarazime na prvek, ktery nalezi spravnemu offsetu
		3 - pokud dorazime na konec offsetu, pak opetovne zatrizovani zaciname na zacatku dalsiho offsetu
	*/
	//alfa pass
	for (u_int i = 0; i < SA_ALFA_LEN; i++) {
		//off[i] to off[i] + f[i]		
		for (u_int j = off[i]; j < off_end[i];j++) {
			//zjisteme symbol na pozici i
			char c = baseMap[m_sequence[s_array[j] + level]];
			u_int c_pos = s_array[j];
			while (c != i) {
				//take the next symbol at particular offset
				u_int sa_entry = s_array[off[c]];
				//store c at correct position
				s_array[off[c]] = c_pos;
				//move offset to the next unsorted position
				off[c]++;
				//reset c
				c = baseMap[m_sequence[sa_entry + level]];
				c_pos = sa_entry;
			}
			s_array[off[i]] = c_pos;
			off[i]++;			
		}
	}
	
	//recursive radix => we don't need to radix 0
	u_int new_start = start + f[0];
	for (u_int i = 1; i < SA_ALFA_LEN; i++) {
		//we also don't need to radix interval of length 1 => it is already sorted
		if (f[i] > 1)
			RadixSortRecursive(new_start, new_start+f[i]-1, level + 1);
		else if (f[i] == 1) {
			sorted++;
			if (sorted % 100000 == 0)
				cout << sorted << endl;
		}

		new_start += f[i];
	}
}

void SuffixArray::SortArray(u_int start, u_int end)
{
	//prepare frequency array
	u_int f[SA_ALFA_LEN];
	for (u_int i = 0; i < SA_ALFA_LEN; i++)
		f[i] = 0;

	//count frequencies in interval
	for (u_int i = start; i <= end; i++)
		f[baseMap[m_sequence[s_array[i]]]]++;

	//offset array - pointers to the beginning of the sa interval for particular symbol
	u_int off[SA_ALFA_LEN];
	u_int off_end[SA_ALFA_LEN];
	off[0] = start;
	off_end[0] = start + f[0];

	u_int tot = start + f[0];
	for (u_int i = 1; i < SA_ALFA_LEN; i++) {
		off[i] = tot;
		off_end[i] = off[i] + f[i];
		tot += f[i];
	}

	/*
	Postup:
	1 - prochazime po jednotlivych offsetech
	2 - umistujeme prvky na spravnou pozici, dokud nenarazime na prvek, ktery nalezi spravnemu offsetu
	3 - pokud dorazime na konec offsetu, pak opetovne zatrizovani zaciname na zacatku dalsiho offsetu
	*/
	//alfa pass
	for (u_int i = 0; i < SA_ALFA_LEN; i++) {
		//off[i] to off[i] + f[i]		
		for (u_int j = off[i]; j < off_end[i]; j++) {
			//zjisteme symbol na pozici i
			char c = baseMap[m_sequence[s_array[j]]];
			u_int c_pos = s_array[j];
			while (c != i) {
				//take the next symbol at particular offset
				u_int sa_entry = s_array[off[c]];
				//store c at correct position
				s_array[off[c]] = c_pos;
				//move offset to the next unsorted position
				off[c]++;
				//reset c
				c = baseMap[m_sequence[sa_entry]];
				c_pos = sa_entry;
			}
			s_array[off[i]] = c_pos;
			off[i]++;
		}
	}

	//recursive radix => we don't need to radix 0
	u_int new_start = start + f[0];
	for (u_int i = 1; i < SA_ALFA_LEN; i++) {
		//omit N
		if(i != 4)
			QuickSort(new_start, new_start + f[i] - 1);
		
		cout << i << "/" << SA_ALFA_LEN << "\t symbols of alphabet already sorted." << endl;
		new_start += f[i];
	}
}

bool SuffixArray::LT(int a, int b)
{
	if (a == b)
		return false;

	u_int id = 0;
	while (true) {
		if (m_sequence[s_array[a] + id] < m_sequence[s_array[b] + id])
			return true;
		else if (m_sequence[s_array[a] + id] > m_sequence[s_array[b] + id])
			return false;
		id++;
	}
	return false;
}

bool SuffixArray::GT(int a, int b)
{
	if (a == b)
		return false;

	u_int id = 0;
	while (true) {
		if (m_sequence[s_array[a] + id] > m_sequence[s_array[b] + id])
			return true;
		else if (m_sequence[s_array[a] + id] < m_sequence[s_array[b] + id])
			return false;
		id++;
	}
	return false;
}

void SuffixArray::QuickSort(u_int left, u_int right) {
	if (left == right)
		return;

	u_int l = left;
	u_int r = right;

	u_int pivot = (l + r) / 2;
	do
	{
		while (LT(l, pivot)) l++;
		while (LT(pivot, r)) r--;
		if (l <= r) {
			u_int aux = s_array[l];
			s_array[l] = s_array[r];
			s_array[r] = aux;

			if (pivot == l)
				pivot = r;
			else if (pivot == r)
				pivot = l;

			l++;
			if (r != 0)
				r--;
		}
	} while (l <= r);

	if (left < r)
		QuickSort(left, r);
	if (l < right)
		QuickSort(l, right);
};

void SuffixArray::Show()
{
	cout << "SA: " << endl;
	for (u_int i = 0; i < m_sequence_length; i++) {
		cout << i << " - " << s_array[i] << " - " << m_sequence[s_array[i]] << endl;
	}
	cout << endl;
}

bool SuffixArray::Verify()
{
	for (u_int i = 1; i < m_sequence_length; i++) {
		//compare alphabetically that s_array[i-1] < s_array[i]
		bool lt = true;
		u_int j = 0;
		while (s_array[i - 1] + j != m_sequence_length && s_array[i] + j != m_sequence_length){
			if (m_sequence[s_array[i - 1]] == 'N' || m_sequence[s_array[i]] == 'N')
				break;

			if (m_sequence[s_array[i - 1] + j] > m_sequence[s_array[i] + j]) {
				lt = false;
				break;
			}
			else if (m_sequence[s_array[i - 1] + j] < m_sequence[s_array[i] + j])
				break;

			j++;			
		}
		if(!lt){
			cerr << "Suffix array incorectly sorted." << endl;
			for (u_int k = s_array[i-1]; k < s_array[i - 1] + 100; k++)
				cout << m_sequence[k];
			cout << endl;
			for (u_int k = s_array[i]; k < s_array[i] + 100; k++)
				cout << m_sequence[k];
			cout << endl;
			Show();
			exit(EXIT_FAILURE);
		}
	}
}

void SuffixArray::Serialize(char* sa_path)
{
	ofstream output(sa_path, ios::binary | ios::app);
	output.write((const char*)s_array, sizeof(u_int)*m_sequence_length);
	output.close();
}

void SuffixArray::SerializeBWT(char* bwt_path)
{
	char* bwt_buf = new char[m_sequence_length];

	for (u_int i = 0; i < m_sequence_length; i++) {
		if (s_array[i] == 0)
			bwt_buf[i] = 0;
		else
			bwt_buf[i] = m_sequence[s_array[i] - 1];
	}

	ofstream output_bwt(bwt_path, ios::binary | ios::app);
	output_bwt.write(bwt_buf, m_sequence_length);
	output_bwt.close();

	delete[] bwt_buf;
}

void SuffixArray::Backtrack(char* sequence, u_int seq_id, u_int& start, u_int& end)
{
	//localize smallest next sym in l column	
	u_int new_start = 1;
	u_int new_end = 0;
	for (u_int i = start; i <= end; i++) {
		if (sequence[seq_id] == m_bwt[i]) {
			new_start = m_offsets[sequence[seq_id]] + m_rank[i];
			break;
		}
	}

	for (u_int i = end; i >= start; i--) {
		if (sequence[seq_id] == m_bwt[i]) {
			new_end = m_offsets[sequence[seq_id]] + m_rank[i];
			break;
		}
	}

	start = new_start;
	end = new_end;

	if (seq_id == 0 || start > end)
		return;
	else
		Backtrack(sequence, --seq_id, start, end);
}

u_int* SuffixArray::Localize(char* src, u_int src_len, u_int& occurences)
{
	u_int start = m_offsets[src[src_len-1]];
	u_int end = m_offsets[src[src_len - 1] + 1] - 1;

	Backtrack(src, src_len - 2, start, end);

	if (start > end) {
		//none found
		occurences = 0;
		return NULL;
	}
	else {
		occurences = end - start + 1;
		u_int* positions = new u_int[occurences];
		for (u_int i = 0; i < occurences; i++) {
			positions[i] = s_array[start + i];
		}

		return positions;
	}
}
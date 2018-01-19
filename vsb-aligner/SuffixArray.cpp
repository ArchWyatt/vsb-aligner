#include "SuffixArray.h"

SuffixArray::SuffixArray(char* sequence, u_int sequence_length)
{
	m_sequence = sequence;
	m_sequence_length = sequence_length;

	s_array = new u_int[sequence_length];
	for (u_int i = 0; i < sequence_length; i++)
		s_array[i] = i;

	RadixSort(0, m_sequence_length - 1, 0);
}

SuffixArray::~SuffixArray()
{
	delete[] s_array;
}

void SuffixArray::RadixSort(u_int start, u_int end, u_int level)
{
	//prepare frequency array
	u_int f[SA_ALFA_LEN];
	for (u_int i = 0; i < SA_ALFA_LEN; i++)
		f[i] = 0;

	//count frequencies in interval
	for (u_int i = start; i <= end; i++)
		f[Utils::BaseMap(m_sequence[s_array[i]+level])]++;

	//offset array - pointers to the beginning of the sa interval for particular symbol
	u_int off[SA_ALFA_LEN];
	off[0] = 0;

	u_int tot = f[0];
	for (u_int i = 1; i < SA_ALFA_LEN; i++) {
		off[i] = tot;
		tot += f[i];
	}

	//make a copy of sa interval
	u_int* sa = new u_int[m_sequence_length];
	memcpy(sa, s_array, sizeof(u_int)*m_sequence_length);

	//radix
	for (u_int i = start; i <= end; i++) {
		//find the proper sa index
		char c = Utils::BaseMap(m_sequence[sa[i] + level]);
		u_int s_index = start + off[c];
		s_array[s_index] = sa[i];		
		off[c]++;
	}
	delete[] sa;

	//recursive radix => we don't need to radix 0
	u_int new_start = start + f[0];
	for (u_int i = 1; i < SA_ALFA_LEN; i++) {
		//we also don't need to radix interval of length 1 => it is already sorted
		if (f[i] > 1)
			RadixSort(new_start, new_start+f[i]-1, level + 1);

		new_start += f[i];
	}
}

void SuffixArray::Show()
{
	cout << "SA: " << endl;
	for (u_int i = 0; i < m_sequence_length; i++) {
		cout << i << " - " << s_array[i] << endl;
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
			for (u_int k = s_array[i-1]; k < m_sequence_length; k++)
				cout << m_sequence[k];
			cout << endl;
			for (u_int k = s_array[i]; k < m_sequence_length; k++)
				cout << m_sequence[k];
			cout << endl;
			exit(EXIT_FAILURE);
		}
	}
}
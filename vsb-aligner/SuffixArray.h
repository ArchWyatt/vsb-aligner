#pragma once

#include "Definitions.h"
#include "Utils.h"

class SuffixArray
{
private:
	/*
		Suffix array container.
	*/
	u_int* s_array;

	/*
		Seuquence to sort.
	*/
	char* m_sequence;
	u_int m_sequence_length;
	
	/*
		Sort suffix array: method recursive radix
	*/
	void RadixSort(u_int start, u_int end, u_int level);

public:
	SuffixArray(char* sequence, u_int sequence_length);
	~SuffixArray();

	//void Serialize(char* sa_path);

	/*
		Checks whether suffix array already exists.
	*/
	//static bool SuffixIndexExists(char* sa_path)

	void Show();

	bool Verify();
};

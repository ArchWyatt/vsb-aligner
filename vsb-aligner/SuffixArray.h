#pragma once

#include <fstream>

#include "Definitions.h"
#include "Utils.h"

using namespace std;

class SuffixArray
{
private:
	/*
		counter of already sorted items.
	*/
	u_int sorted = 0;

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
		Alfa character map.
	*/
	char* baseMap;
	
	/*
		Sort suffix array: method recursive radix
	*/
	void RadixSortRecursive(u_int start, u_int end, u_int level);

	/*
		First split the entries of an array by radix into groups beginning with the same symbol. Each symbol other than N then quicksort.
	*/
	void SortArray(u_int start, u_int end);

	/*
		Quicksort suffix array.
	*/
	void QuickSort(u_int left, u_int right);
	/*
		Compare string at two suffix indexec a and b.
		returns true if a < b
	*/
	bool LT(int a, int b);
	/*
		returns true if a > b
	*/
	bool GT(int a, int b);

public:
	SuffixArray(char* sequence, u_int sequence_length);
	~SuffixArray();

	/*
		Will serialize sa.
	*/
	void Serialize(char* sa_path);

	/*
		Will serialize as bwt.
	*/
	void SerializeBWT(char* bwt_path);
	
	void Show();

	bool Verify();
};

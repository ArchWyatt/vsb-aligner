#pragma once

#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <string>

#include "Definitions.h"
#include "GenomeRegion.h"
//#include "Interval.h"

using namespace std;

static const u_int INTERVAL_ARGUMENTS = 6;
static const u_int INTERVAL_ARGUMENTS_GENE_FUSION = 7;
static const u_int GENOME_REGION_ARGUMENTS = 5;

class Utils {
public:
	static char* ParseArgByDoubleSpace(char* src, u_int src_len, u_int id, u_int& arg_len) {
		u_int item_idx = 0;
		u_int ptr = 0;
		while (item_idx != id){
			while (src[ptr++] != ' ') {
				if (ptr == src_len) {
					//not found
					return NULL;
				}
			}
			//move one ahead i.e. over the next space
			ptr++;
			item_idx++;
		}

		u_int aux = ptr;
		while (src[aux] != ' ') aux++;
		//ptr aux points to next space
		arg_len = aux - ptr;

		char* c_word = new char[arg_len + 1];
		c_word[arg_len] = 0;
		memcpy(c_word, src + ptr, arg_len);

		return c_word;
	};

	/* tag:value
		- returns value in u_int
	*/
	static u_int TagParsePosition(char* tag, u_int tag_len) {
		u_int ptr = 0;
		while (ptr < tag_len && tag[ptr++] != ':') {}

		//tag[ptr] contains first symbol of value
		//tag_len is one above than strlen
		//u_int value_len = tag_len - ptr;

		return atol(tag + ptr);
	};


	static GenomeRegion* LineToRegion(char* line)
	{
		u_int line_len = strlen(line);
		if (line_len == 0)
			return NULL;

		char* args[GENOME_REGION_ARGUMENTS];
		u_int c = 0;
		u_int i = 0;
		while (line[i] != 0 && c < GENOME_REGION_ARGUMENTS) {
			//scan for word delimited by tab			
			u_int l = i;
			while (!(line[l] == '\t' || line_len == l)) l++;
			//l => length of the word
			char* word = new char[l - i + 1];
			for (u_int j = 0; j < l - i; j++) {
				word[j] = line[i + j];
			}
			word[l - i] = 0;
			i = l + 1;
			args[c++] = word;
		}

		if (c != GENOME_REGION_ARGUMENTS) {
			cerr << "ERROR: reading not enough arguments for region." << endl;
			exit(EXIT_FAILURE);
		}

		GenomeRegion* gen_region = new GenomeRegion(args[0], atoi(args[1]), stoul(args[2]), atoi(args[3]), atoi(args[4]));

		for (u_int i = 1; i < GENOME_REGION_ARGUMENTS; i++) delete[] args[i];

		return gen_region;
	};

	/* Does target begins with pattern? */
	static bool BeginsWith(char* pattern, char* target)
	{
		if (strlen(target) < strlen(pattern))
			return false;

		for (u_int i = 0; i < strlen(pattern); i++) {
			if (pattern[i] != target[i]) {
				return false;
			}
		}
		return true;
	};

	/* Does the target ends with pattern? */
	static bool EndsWith(char* pattern, char* target)
	{
		u_int t_len = strlen(target);
		u_int p_len = strlen(pattern);
		if (t_len < p_len)
			return false;

		for (u_int i = 1; i <= p_len; i++) {
			if (target[t_len - i] != pattern[p_len - i])
				return false;
		}
		return true;
	};

	/* String copy - expects c-style string*/
	static char* StrCopy(char* src)
	{
		u_int s_len = strlen(src);
		char* new_s = new char[s_len + 1];
		new_s[s_len] = 0;
		memcpy(new_s, src, sizeof(char)*s_len);
		return new_s;
	};

	/* String append */
	static char* StrAppend(char* src, const char* app_src)
	{
		u_int s_len = strlen(src) + strlen(app_src);
		char* new_s = new char[s_len + 1];
		new_s[s_len] = 0;
		memcpy(new_s, src, sizeof(char)*strlen(src));
		memcpy(new_s + strlen(src), app_src, sizeof(char)*strlen(app_src));

		return new_s;
	};

	/*
		String uppercase
		TODO: is it faster when checked for <a,z> and add (a-A)?
	*/
	static void StrUppercase(char* src, u_int s_len)
	{
		for (u_int i = 0; i < s_len; i++) {
			if (src[i] == 'a')
				src[i] = 'A';
			else if (src[i] == 'c')
				src[i] = 'C';
			else if (src[i] == 'g')
				src[i] = 'G';
			else if (src[i] == 't')
				src[i] = 'T';
			else if (src[i] == 'n')
				src[i] = 'N';
		}
	};

	/*
		String reverse - returns new string
	*/
	static char* StrReverse(char* src)
	{
		u_int s_len = strlen(src);
		char* new_s = new char[s_len + 1];
		new_s[s_len] = 0;

		for (u_int i = 0; i < s_len; i++)
			new_s[s_len - i - 1] = src[i];
		
		return new_s;
	};

	/*
		Function that matches pattern in the position into target splitted by delimiter.
	*/
	static bool IsInPositionByDelimiter(char* pattern, char* target, char delimiter, u_int position)
	{
		u_int t_len = strlen(target);
		u_int p_len = strlen(pattern);

		if (t_len < p_len)
			return false;

		u_int t_ptr = 0;
		u_int delim_no = 0;
		while (delim_no != position && t_ptr < t_len) {
			if (target[t_ptr] == delimiter)
				delim_no++;

			t_ptr++;
		}

		if (t_ptr == t_len && p_len != 0)
			return false;
		else if(t_ptr == t_len && delim_no == position)
			return true;

		for (u_int i = 0; i < p_len; i++) {
			if (target[t_ptr + i] != pattern[i])
				return false;
		}
		return true;
	};

	static char* Substring(char* target, u_int from, u_int length)
	{
		char* new_s = new char[length + 1];
		new_s[length] = 0;

		memcpy(new_s, target + from, sizeof(char)*length);

		return new_s;
	};

	/*
		Maps bases into the range 0-5
	*/
	static u_char BaseMap(char base) {
		if (base == 0)
			return 0;
		else if (base == 'A')
			return 1;
		else if (base == 'C')
			return 2;
		else if (base == 'G')
			return 3;
		else if (base == 'N')
			return 4;
		else if (base == 'T')
			return 5;
		else {
			cerr << "Unknow base: " << base << endl;
			exit(EXIT_FAILURE);
		}
	}
};

#endif
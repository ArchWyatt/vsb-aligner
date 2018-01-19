#pragma once

/*
	Hash
		- default implementation of hash - each item of hash array is assigned with a list
		- key => char*
*/

#include <math.h>
#include <cstring>

#include "Definitions.h"
#include "List.h"

static const u_int DEFAULT_HASH_CAPACITY = 256;

using namespace std;

template<class T>
class HashItem
{
public:
	HashItem(char* new_key, T* new_value) {
		key = new_key;
		value = new_value;
	};
	~HashItem() {
		delete[] key;
	};

	char* key;
	T* value;
};

template<class T>
class Hash
{
private:
	u_int capacity;
	/*
		array to hold hash objects
	*/
	List<HashItem<T>>** h_array;

public:
	Hash(u_int new_capacity) {
		capacity = new_capacity;

		h_array = new List<HashItem<T>>*[capacity];
		//initialize each item with empty list
		for (u_int i = 0; i < capacity; i++)
			h_array[i] = new List<HashItem<T>>();		
	};

	/*
		Default constructor with default Hash array capacity.
	*/
	Hash() : Hash(DEFAULT_HASH_CAPACITY){};	

	~Hash() {
		delete[] h_array;
	};

	/*
		Bracket operator overloading.
	*/
	// returns value
	T* Get(char* key){
		//compute hash value
		u_int hval = HashFunction(key);
		//access proper list indexed by key
		List<HashItem<T>>* proper_list = h_array[hval];
		//check if the proper is list is not empty
		if (proper_list->Length() == 0) return NULL;
		//iterate over the list and find the proper object
		ListIterator<HashItem<T>> proper_list_iterator(proper_list->First());

		while (proper_list_iterator.Current() != NULL) {
			//vytahneme klic z iteratoru
			char* proper_key = proper_list_iterator.Current()->Value()->key;
			if (strcmp(proper_key, key) == 0) {
				return proper_list_iterator.Current()->Value()->value;
			}

			proper_list_iterator.Next();
		}

		//if we are here then there is no value corresponding to the key
		return NULL;
	};

	/*
		Add values to the hash
	*/
	void Add(char* key, T* value)
	{
		//compute hash value
		u_int hval = HashFunction(key);
		//access proper list indexed by key
		List<HashItem<T>>* proper_list = h_array[hval];
		//create new hash item
		HashItem<T>* new_item = new HashItem<T>(key, value);
		proper_list->Append(new_item);
	}

	/*
		Simple hash function
	*/
	u_int HashFunction(char* key) {
		u_int t = 1;

		for (u_int i = 0; i < strlen(key); i++) {
			t += powl(key[0], i+2);
		}

		return t % capacity;
	}
};


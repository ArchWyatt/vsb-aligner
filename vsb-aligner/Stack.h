#pragma once

#include "Definitions.h"

struct StackItem
{
	u_int start;
	u_int end;
	u_int level;

	StackItem* prev;
	StackItem* next;	
};

class Stack
{
private:
	StackItem* root;

public:
	Stack() {
		root = NULL;
	};

	~Stack() {
		delete root;
	};

	void Push(u_int start, u_int end, u_int level) {
		StackItem* item = new StackItem();
		item->start = start;
		item->end = end;
		item->level = level;
	};

	StackItem* Pop() {

	};
};

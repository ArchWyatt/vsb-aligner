#pragma once
#ifndef __LIST_H__
#define __LIST_H__

template<class T> class ListItem
{
private:
	T* value;
	ListItem* next;
	ListItem* prev;

public:
	ListItem(T* arg_val) {
		value = arg_val;
		next = NULL;
		prev = NULL;
	};

	~ListItem() {};

	T* Value() {
		return value;
	};

	ListItem<T>* Next() {
		return next;
	};

	ListItem<T>* SetNext(ListItem<T>* new_next) {
		next = new_next;
		return next;
	};

	bool HasNext() {
		if (next == NULL)
			return false;

		return true;
	};

	ListItem<T>* Prev() {
		return prev;
	};

	ListItem<T>* SetPrev(ListItem<T>* new_prev) {
		prev = new_prev;
		return prev;
	};

	bool HasPrev() {
		if (prev == NULL)
			return false;

		return true;
	};
};

template<class T> class ListIterator
{
private:
	ListItem<T>* current;
	ListItem<T>* first;	
	ListItem<T>* last_used;

public:
	ListIterator() {
		current = NULL;
		last_used = NULL;
	};

	ListIterator(ListItem<T>* current_item) {
		first = current_item;
		current = current_item;
	};

	~ListIterator() {};

	ListItem<T>* Current() {
		return current;
	};

	void SetCurrent(ListItem<T>* new_current) {
		current = new_current;
		last_used = NULL;
	};

	void Next() {
		last_used = current;
		current = current->Next();
	};

	bool HasNext() {
		return current->HasNext();
	};

	void Prev() {
		last_used = current;
		current = current->Prev();
	};

	bool HasPrev() {
		return current->HasPrev();
	};

	ListItem<T>* LastUsed() {
		return last_used;
	};

	void Reset() {
		current = first;
	};
};

template<class T> class List
{
private:
	ListItem<T>* root;
	ListItem<T>* last;

	u_int length;

public:
	List() {
		root = NULL;
		last = NULL;
		length = 0;
	};

	~List() {
		ListIterator<T> list_iterator(root);

		while (list_iterator.Current() != NULL) {
			ListItem<T>* current = list_iterator.Current();
			list_iterator.Next();
			delete current;
		}
	};

	u_int Length() {
		return length;
	};

	ListItem<T>* First() {
		return root;
	};

	ListItem<T>* Last() {
		return last;
	};

	void Append(T* value) {
		if (root == NULL){
			last = new ListItem<T>(value);
			root = last;
		}
		else {
			ListItem<T>* aux = last;
			last = aux->SetNext(new ListItem<T>(value));
			last->SetPrev(aux);
		}
		length++;
	};

	void Insert(ListItem<T>* prev, ListItem<T>* cur) {
		if (prev == NULL) {
			//place it as a root
			ListItem<T>* aux = root;
			root = cur;
			root->SetNext(aux);
		}
		else {
			ListItem<T>* next_aux = prev->Next();
			prev->SetNext(cur);
			cur->SetNext(next_aux);
		}
	};

	T* Find(T* value) {
		if (root == NULL)
			return NULL;

		ListItem<T>* aux = root;
		while (aux != NULL) {
			if (aux->Value() == value) {
				return aux->Value();
			}
			aux = aux->Next();
		}
		return NULL;
	};

	/* return true if the value is removed otherwise false */
	bool Remove(T* value) {
		if (root == NULL)
			return false;

		ListItem<T>* aux = root;
		ListItem<T>* aux_last = NULL;
		while (aux != NULL) {
			if (aux->Value() == value) {
				if (aux_last != NULL) {
					aux_last->SetNext(aux->Next());					
				}

				if (value == root->Value())
					root = aux->Next();

				delete aux;
				length--;
				return true;
			}
			aux_last = aux;
			aux = aux->Next();
		}
		return false;
	};

	/* item is removed based on the knowledge of the last item */
	bool Remove(T* value, ListItem<T>* last_item){
		if (last_item == NULL) {
			ListItem<T>* aux = root;
			if (aux->Value() != value)
				return false;

			if (root == last)
				last = NULL;
			
			root = aux->Next();
			delete aux;
			length--;
			return true;
		}
		else {
			ListItem<T>* aux = last_item->Next();
			if (aux->Value() != value)
				return false;

			last_item->SetNext(last_item->Next()->Next());
			if (last_item->Next() == NULL)
				last = last_item;

			delete aux;
			length--;
			return true;
		}
		return false;
	};
};

#endif
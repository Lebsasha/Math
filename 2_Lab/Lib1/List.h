#ifndef LIST_H
#define LIST_H
#include "stdafx.h"
#include <iostream>
#include <string>
using namespace std;
struct Student
{
	string Name;
	string Surname;
	string SecName;
	int Ages;
	int Male;
	int Kurs;
	int Marks;
	void View (void)
	{
		cout<<this->Surname;
		cout<<endl<<this->Name;
		cout<<endl<<this->SecName;
		cout<<endl<<this->Ages;
		cout<<endl<<this->Male;
		cout<<endl<<this->Kurs;
		cout<<endl<<this->Marks;
		cout<<endl;
	}
	bool Compare (Student& Somewho)
	{
		if (this->Surname.compare (Somewho.Surname) >= 0)
			return 1;
		return 0;
	}
};
struct Node
{
	Node* pNext;
	Student Pal;
	Node (const Student& Temp): pNext(NULL), Pal (Temp)
	{}
	Node (): pNext(NULL)
	{}
};
bool Read_from_t_file_list_st (ifstream& iFile, Node*& Head);
void View (const Node* Head);
void Delete (Node* Prev);
void Add (Node* Head, Node* Some);
void Remove (Node* Prev);
Node* Goto (Node* Head, int co);



struct Node1
{
	Node1* pNext;
	Node1* pPrev;
	double Pal;
	Node1 (): Pal(0), pNext(NULL), pPrev(NULL)
	{}
	Node1 (double x): Pal(x), pNext(NULL), pPrev (NULL)
	{}
};
void Add (Node1* Place, Node1* Some);
void Delete (Node1* Prev);
double Calc (Node1* Head, Node1* Tail);
#endif

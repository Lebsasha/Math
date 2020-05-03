#ifndef LIST_H
#define LIST_H
#include <iostream>
#include <string>
#if __cplusplus < 201103L
#define nullptr NULL
#endif
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
    Student (void): Name(), Surname(), SecName(), Ages(0), Male(0), Kurs(0), Marks(0)
    {}
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
    Node (const Student& Temp): pNext(nullptr), Pal (Temp)
    {}
    Node (): pNext(nullptr), Pal()
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
    Node1 (): pNext(nullptr), pPrev(nullptr), Pal(0)
    {}
    Node1 (double x): pNext(nullptr), pPrev (nullptr), Pal(x)
    {}
};
void Add (Node1* Place, Node1* Some);
void Delete (Node1* Prev);
double Calc (Node1* Head, Node1* Tail);
#endif

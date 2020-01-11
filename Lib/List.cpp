#include "List.h"
#include <iostream>
#include <fstream>
using namespace std;
void View (const Node* Head)
{
    for (Node* p = Head->pNext; p; p = p->pNext)
    {
        p->Pal.View();
        cout<<endl;
    }
    return;
}
void Read_from_f_txt (ifstream& iFile, Student* Stu)
{
    iFile>>Stu->Surname;
    iFile>>Stu->Name;
    iFile>>Stu->SecName;
    iFile>>Stu->Ages;
    iFile>>Stu->Male;
    iFile>>Stu->Kurs;
    iFile>>Stu->Marks;
    return;
}
void Add (Node* Head, Node* Some)
{
    for (Node* p = Head; p; p = p->pNext)
    {
        if (!p->pNext || p->pNext->Pal.Compare (Some->Pal))
        {
            Some->pNext = p->pNext;
            p->pNext = Some;
            break;
        }
    }
    return;
}
void Delete (Node* Prev)
{
    Node* Temp = Prev->pNext;
    Prev->pNext = Prev->pNext->pNext;
    delete Temp;
    return;
}
void Remove (Node* Prev)
{
    Prev->pNext = Prev->pNext->pNext;
    return;
}
bool Read_from_t_file_list_st (ifstream& iFile, Node*& Head)
{
<<<<<<< HEAD
    Node* One = new Node;
    Read_from_f_txt(iFile, &One->Pal);
    Head->pNext = One;
//Node* Prev = Head; // Unused
    while (!iFile.fail())
    {
        One = new Node;
        if (One)
        {
            Read_from_f_txt(iFile, &One->Pal);
            Add (Head, One);
        }
        else
        {
            return 0;
        }
    }
    return 1;
=======
Node* One = new Node;
Read_from_f_txt(iFile, &One->Pal);
Head->pNext = One;
Node* Prev = Head;
while (!iFile.fail())
{
One = new Node;
if (One)
{
Read_from_f_txt(iFile, &One->Pal);
Add (Head, One);
}
else
{
return 0;
}
}
return 1;
>>>>>>> ec95f91ed173f170c4c98998802d56dbf22b6d2e
}
Node* Goto (Node* Head, int co)
{
    while (co--)
    {
        Head = Head->pNext;
    }
    return Head;
}
/*bool Move (Node* Prev_To, Node* Prev_From)
{
Prev_To->pNext = Prev_From->pNext;
Prev_From->pNext = Prev_To->pNext->pNext;
Prev_To->pNext->pNext =
}*/
void Add (Node1* Place, Node1* Some)
{
    Place->pNext = Some;
    Some->pPrev = Place;
    return;
}
void Delete (Node1* Prev)
{
    Node1* Temp = Prev->pNext;
    Prev->pNext = Prev->pNext->pNext;
    Prev->pNext->pPrev = Prev;
    delete Temp;
    return;
}
double Calc (Node1* Head, Node1* Tail)
{
<<<<<<< HEAD
    double Res = 0;
    const Node1* Head1 = Head;
    while ( Head1 != Tail)
    {
        Res *= Head->Pal+Tail->Pal;
        Head = Head->pNext;
        Tail = Tail->pPrev;
    }
    return Res;
}
=======
	double Res = 0;
	const Node1* Head1 = Head;
	while ( Head1 != Tail)
	{
		Res *= Head->Pal+Tail->Pal;
		Head = Head->pNext;
		Tail = Tail->pPrev;
	}
	return Res;
}
>>>>>>> ec95f91ed173f170c4c98998802d56dbf22b6d2e

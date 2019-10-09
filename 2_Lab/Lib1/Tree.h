#include <iostream>
#include <time.h>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iterator>
template <class T>
struct Nt
{
    T Pal;
    Nt* pLeft;
    Nt* pRight;
    Nt(void): Pal(), pLeft(NULL), pRight(NULL)
    {}
    //Nt (T& a): Pal(a), pLeft(NULL), pRight(NULL)
    //{}
    Nt (T a): Pal(a), pLeft(NULL), pRight(NULL)
    {}
    ~Nt (void)
    {Pal = 0;}
        //delete Pal;
};
template<class T>
class Binary_tree
{
    Nt<T>* Root;
    Nt<T>* Temp;
    int N_le1 (Nt<T>*& Temp1) const
    {
        int Res = 0;
        if (Temp1 == NULL)
            return 0;
        if (Temp1->pLeft == NULL && Temp1->pRight == NULL)
            return 1;
        else
            Res = N_le1(Temp1->pLeft) + N_le1(Temp1->pRight);
        return Res;
    }
    void Add1 (Nt<T>*& t, const T& a)
    {
        if (!t)
        {
            t = new Nt<T>(a); ++i;
            if (t)
                return;
            return;
        }
        else
        {
            //Nt<T>* prevTemp = Temp;
            if (a < t->Pal)
                Add1 (t->pLeft, a);
            else
                Add1 (t->pRight, a);
        }
    }
    void View1 (Nt<T>* Temp1, const bool b) const
    {
        if (Temp1)
        {
            View1 ((b) ? Temp1->pLeft : Temp1->pRight, b);
            cout<<Temp1->Pal<<' ';
            View1 ((b) ? Temp1->pRight : Temp1->pLeft, b);
        }
        return;
    }
    void Copy (Nt<T>*& pTree, Nt<T>*& pBTree)
    {
        if (!pBTree)
        {
            return;
        }
        pTree = new Nt<T> (pBTree->Pal); ++i;
        if (pTree)
        {
            Copy (pTree->pLeft, pBTree->pLeft);
            Copy (pTree->pRight, pBTree->pRight);
        }
    }
public:
    Binary_tree (T& a): Root (new Nt<T> (a)), Temp (Root)
    {++i;}
    Binary_tree (void): Root(new Nt<T>(0)), Temp(Root)
    {++i;}
    Binary_tree (Binary_tree& Tr): Root (NULL), Temp(Root)
    {
        //Nt<T>* pA = Tr.Root;
        Copy (Root, Tr.Root);
    }
    bool Add (const T& a)
    {
        Add1 (Root, a);
        return 1;
    }
    void Delete (const T& a)
    {
        Nt<T>* Temp1 = Temp;
        while (Temp =  Find_Prev (a))
        {
            if (Temp->pLeft->Pal == a)
            {
                Temp1 = Temp->pLeft;
                Temp->pLeft = 0;
            }
            else
            {
                Temp1 = Temp->pRight;
                Temp->pRight = 0;
            }
            delete Temp1; --i;
        }
        return;
    }
    Nt<T>*& Find_Prev (const T& a) const
    {
        Nt<T>* PrevT = 0;
        Nt<T>* Temp1 = Root;
        while ( Temp1 && Temp1->Pal != a)
        {
            PrevT = Temp1;
            if (a < Temp1->Pal)
                Temp1 = Temp1->pLeft;
            else
                Temp1 = Temp1->pRight;
        }
        if (Temp1)
            return PrevT;
        return Temp1;
    }
    int N_le(void) const
    {
        int Res = N_le1(Root->pLeft) + N_le1(Root->pRight);
        return Res;
    }
    void Del_all_aft (Nt<T>*& Tem)
    {
        if (Tem == NULL)
            return;
        Del_all_aft(Tem->pLeft);
        Del_all_aft(Tem->pRight);
        if (Tem->pLeft == NULL && Tem->pRight == NULL)
        {
            delete Tem; --i;
            Tem = NULL;
            return;
        }
        return;
    }
    void View (const bool b = 0) const
    {
        View1 (Root, b);
        return;
    }
    ~Binary_tree (void)
    {
        Del_all_aft(Root);
        //cout<<i<<"\n";
    }
};

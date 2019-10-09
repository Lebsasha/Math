#include "stdafx.h"
//#include "Vector.h"
#include <iostream>
using namespace std;
/*Vector (const Physics** pA, const int N)
{
Physics** pa = pA + N -1;
a = new Physics* [N];
for (--a; pA <= pa; ++pA)
{
*(++a) = *pA;
}
}*/
/*void Vector::Rand (int Biggest_Num = 5)
{
++Biggest_Num;
for (Physics** pA = a + N-1; pA >= a; --pA)
*pA = rand()%Biggest_Num;
}
friend ostream& Vector::operator<< (ostream& a, const Vector& v)
{
Physics** End = v.a + v.N - 1;
for (Physics** pa = v.a; pa <= End;++pa)
a<<*pa;
return a;
}
friend istream& Vector::operator>> (istream& a, Vector& v)
{
	//cin>>v.N;
	//pEnd = v.a + v.N - 1;
	//for (Physics** pa = v.a; pa <= pEnd; ++pa)
	//{
	//	cin>>*pa;
	//}
	return a;
}
inline void Vector::Set_size (int N)
{
this->N = N;
return;
}
inline int Vector::Get_size (void)
{
return N;
}
Physics*& Vector::operator[] (const int i)
{
if (i < N)
return *(a+i);
else
exit (1);
}
Physics* Vector::operator+ (Vector& a)
{
const int N_min = (N < a.N ) N : a.N;
for (Physics** pA = this->a + N_min - 1, **pB = a.a + N_min - 1; pA <= this->a; --pA, --pB)
{
*pA += *pB;
}
return *this;
}*/

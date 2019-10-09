#include <time.h>
#include <iostream>
#include <cstdlib>
#include <iomanip>
using namespace std;
void randoe (int* pArray, const int Nstr, const int Nstb, int biggest_1) /*randoe (&A[0][0], n, n, 50);*/
{
    int i = 0;
    srand(static_cast<int>(time(NULL)));
    while (i++ < Nstb*Nstr)
        {
        *(pArray++) = rand()%biggest_1;
        }
        return;
}
void randoe (double* pArray, const int Nstr, const int Nstb, int biggest_1) /*randoe (&A[0][0], n, n, 50);*/
{
    int i = 0;
    srand(static_cast<int>(time(NULL)));
    while (i++ < Nstb*Nstr)
        {
        *(pArray) = rand()%biggest_1+rand()%100/100.0;
        }
        return;
}
void randoe (float* pArray, const int Nstr, const int Nstb, int biggest_1) /*randoe (&A[0][0], n, n, 50);*/
{
    int i = 0;
    srand(static_cast<int>(time(NULL)));
    while (i++ < Nstb*Nstr)
        {
        *(pArray++) = rand()%biggest_1+rand()%100/100.0;
        }
        return;
}
void View (int* pArray, int Nstr, int Nstb)
{
    while (Nstr--)
        {
        for (int j = 0; j < Nstb-1; j++)
        {
        cout<<setw(3)<<(*pArray)<<" ";
        pArray++;
        }
		cout<<setw(3)<<*pArray++;
        cout<<endl;
        }
    cout<<endl;
    return;
}
void View (double* pArray, int Nstr, const int Nstb)
{
    while (Nstr--)
        {
        for (int j = 0; j < Nstb-1; j++)
        {
        cout<<setw(3)<<(*pArray)<<" ";
        pArray++;
        }
		cout<<setw(3)<<*pArray++;
        cout<<endl;
        }
    cout<<endl;
    return;
}
void View (float* pArray, int Nstr, int Nstb)
{
    while (Nstr--)
        {
        for (int j = 0; j < Nstb-1; j++)
        {
        cout<<setw(3)<<(*pArray)<<" ";
        pArray++;
        }
		cout<<setw(3)<<*pArray++;
        cout<<endl;
        }
    cout<<endl;
    return;
}
void View (double** A, const int n, int m)
	{
	if (!m)
	{
		m = n;
	}
	for (int i = 0; i < n; ++i)
		{
		View (*(A++), 1, m);
		}
		return;
	}
void View (float** A, const int n, int m)
	{
	if (!m)
	{
		m = n;
	}
	for (int i = 0; i < n; ++i)
		{
		View (*(A++), 1, m);
		}
		return;
	}
void View (int** A, const int n, int m)
	{
	if (!m)
	{
		m = n;
	}
	for (int i = 0; i < n; ++i)
		{
		View (*(A++), 1, m);
		}
		return;
	}
#include <vector>
#include <iterator>
void View (vector<double>& El1)
{
    for (vector<double>::const_iterator i = El1.begin(); i < El1.end() - 1; ++i)
        cout<<*i<<' ';
    cout<<*(El1.end()-1);
    return;
}
void Nul(double* p, const int n)
{
	double* pb = p;
	double* pe = p + n;
	while (pb != pe)
	{
	*(pb++) = 0;
	}
	return;
}
void Nul(int* p, const int n)
{
	int* pb = p;
	int* pe = p + n;
	while (pb != pe)
	{
	*(pb++) = 0;
	}
	return;
}
void Nul(float* p, const int n)
{
	float* pb = p;
	float* pe = p + n;
	while (pb != pe)
	{
	*(pb++) = 0;
	}
	return;
}
void Nul(char* p, const int n)
{
	char* pb = p;
	char* pe = p + n;
	while (pb != pe)
	{
	*(pb++) = 0;
	}
	return;
}
bool Sim (int n)
{
	for (int i = 2;i <= n - 2; ++i)
		{
			if (!(n%i))
			{
				return 0;
			}
		}
		return 1;
}

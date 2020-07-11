#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <cstdlib>
/*void Null (int* T, const int N)
{
for (int*t = T+N-1; t >= T; --t)
    *t = 0;
return;
}*/
template <class T>
class Vector
{
    int N;
    T* a;
public:
    Vector (): N(0), a(NULL)
    {}
    explicit Vector (const T& A): N(1), a(new T [N])
    {
        *a = A;
    }
    /*Vector (const T* pA, const int N)
    {
    T* pa = pA + N -1;
    a = new T [N];
    for (--a; pA <= pa; ++pA)
    {
    *(++a) = *pA;
    }
    }*/
    Vector (const Vector& A): N(A.N), a(new T[N])
    {
        --a;
        T* End = A.a + A.N-1;
        for (T* pA = A.a; pA <= End; ++pA)
        {
            *(++a) = *pA;
        }
    }
    void Rand (int Biggest_Num = 5)
    {
        ++Biggest_Num;
        for (T* pA = a + N-1; pA >= a; --pA)
            *pA = rand()%Biggest_Num;
    }
    friend std::ostream& operator<< (std::ostream& str, const Vector& v)
    {
        T* End = v.a + v.N - 1;
        for (T* pa = v.a; pa <= End; ++pa)
        {
            (*pa)->View();///@attention Why?
            str << std::endl;
        }
        return str;
    }
    inline void Set_size (const int M)
    {
        this->N = M;
   }
    inline int Get_size () const
    {
        return N;
    }
    Vector operator+ (const Vector& vA)
    {
        const int N_min = (N < vA.N ) ? N : vA.N;
        for (Physics** pA = this->a + N_min - 1, **pB = vA.a + N_min - 1; pA >= this->a; --pA, --pB)
        {
            *pA += *pB;
        }
        return this;
    }
    bool Add (const T& El)
    {
        T* M = new T[N+1];
        if (!M)
            return false;
        if (!a)
        {
            *M = El;
        }
        else
        {
            for(T* pA = a + N - 1, *pB = M + N; pB >= M; --pA, --pB)
            {
                if (pA >= a && (*pA)->Compare(El))
                    *pB = *pA;
                else
                {
                    *(pB--) = El;
                    for (; pA >= a && pB >= M; --pA, --pB)
                        *pB = *pA;
                }
            }
        }
        delete[] a;
        ++N;
        a = M;
        return true;
    }
    void Delete (const int i)
    {
        if (i > N || i <= 0)
            return;
        T Mem = *(a + i - 1);
        delete Mem;
        T* pEnd = a + --N;
        for (T* pA = a + i - 1; pA < pEnd; ++pA)
            *pA = *(pA + 1);
   }
    int Export (std::ofstream& osa)
    {
        osa.write(reinterpret_cast<char*>(&N), sizeof(int));
        T* pEnd = a + N;
        for (T* pA = a; pA < pEnd; ++pA)
            (*pA)->Export(osa);
        return N;
    }
    int Import (std::ifstream& isa)
    {
        if (isa.eof())
            return 0;
        isa.read(reinterpret_cast<char*>(&N), sizeof(int));
        delete[] a;
        a = new T[N];
        char* C = new char [4];
        T* pEnd = a + N;
        for (T* pA = a; pA < pEnd; ++pA)
        {
            if (isa.eof())
                return 0;
            isa.read(C, 4);
            std::string S(C);
            if (S == "Gra")
                *pA = new Gravity;
            else if (S == "EMI")
                *pA = new EMI;
            else if (S == "Mat")
                *pA = new Matter;
            else if (S == "Str")
                *pA = new Strong;
            else if (S == "Wea")
                *pA = new Weak;
            else
            {
                std::cout<<"Wrong F.bin file"<<std::endl;
                return 0;
            }
            if (!(*pA)->Import(isa))
                return 0;
        }
        delete[] C;
        return N;
    }
    T& operator[] (const int i) const
    {
        if (i < N)
            return *(a+i);
        else
            exit (1);
    }
    virtual ~Vector ()
    {
        T* pEnd = a + N;
        for (T* pA = a; pA < pEnd; ++pA)
            delete *pA;
        delete[] a;
    }
};
/*class Vector
{
Physics** a;
int N;
public:
    Vector ()
    {
    N = 1;
    a = new Physics* [N];
    *a = 0;
    }
    Vector (const Physics* A)
    {
    N = 1;
    a = new Physics* [N];
    *a = A;
    }
    Vector (Physics** pA, const int N) //-??????????????
    {
    Physics** pa = pA + N -1;
    a = new Physics* [N];
    for (--a; pA <= pa; ++pA)
    {
    *(++a) = *pA;
    }
    }
    Vector (Vector& A)
    {
    N = A.N;
    a = new Physics* [N];
    --a;
    Physics** End = A.a + A.N-1;
    for (Physics** pA = A.a; pA <= End; ++pA)
    {
    *(++a) = *pA;
    }
    }

    void Rand (int Biggest_Num = 5);
    friend ostream& operator<< (ostream& a, const Vector& v);
	friend istream& operator>> (ostream& a, Vector& v);
    inline void Set_size (int N);
    Physics* operator+ (Vector& a);
    Physics*& operator[] (const int i);
	~Vector (void)
	{
	delete[] a;
	}
};*/
#endif

#ifdef VISUAL_STUDIO
#include "stdafx.h"
#endif // VISUAL_ST
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cfloat>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "Different.h"
#if __cplusplus < 201103L
#define nullptr NULL
#endif

//TODO Matrix move ctor
template <class T>
class Matrix
{
    Matrix plus_or_minus (const Matrix& addend, bool if_minus)
    {
            assert (N == addend.N && M == addend.M );
            T* p_ans = new T[N * M];
            T* p_end = addend.p_data + addend.N * addend.M - 1;
            if (!if_minus)
                for (T* p_a = p_data - 1, *p_b = addend.p_data - 1, *p_c = p_ans - 1; p_b < p_end;)
                {
                    *(++p_c) = *(++p_a) + *(++p_b);
                }
            else
                for (T* p_a = p_data - 1, *p_b = addend.p_data - 1, *p_c = p_ans - 1; p_b < p_end;)
                {
                    *(++p_c) = *(++p_a) - *(++p_b);
                }
            Matrix answer (p_ans, N, M);
            delete[] p_ans;
            return answer;
    }
    void increase_name_count ()
    {
        std::string a;
        for (int i = ++count; i > 1; i /= 10)
        {
            a += '0' + i%10;
        }
        reverse (a.begin(), a.end());
        name += a;
    }
protected:
    std::string name;
    int N;
    int M;
    T* p_data;
    static int count;
public:
    explicit Matrix (bool if_null = true): name("Matrix "), N(1), M(1), p_data (new T [N * M])
    {
        if (if_null)
            fill_nulls();
        increase_name_count();
    }
    Matrix (T* pointer, const int a, const int b): name("Matrix "), N(a), M(b), p_data(new T [N * M])
    {
        if (pointer)
        {
            for (T* p_ptr_curr = pointer + N * M - 1,* p_curr = p_data + N * M-1; p_ptr_curr >= pointer; --p_ptr_curr,--p_curr)
            {
                *(p_curr) = *p_ptr_curr;
            }
        }
        increase_name_count();
    }
    Matrix (const int a, const int b, bool if_null = true): name("Matrix "), N(a >= 0 ? a : 0), M(b >= 0 ? b : 0), p_data(new T [N * M])
    {
        if (if_null)
            fill_nulls();
        increase_name_count();
    }
    Matrix (const Matrix& matrix): name("Matrix "), N(matrix.N), M(matrix.M), p_data(nullptr)
    {
        *this = matrix;
        increase_name_count();
    }
//    Matrix(Matrix&& mmty) noexcept: name(std::move(mmty.name)), N(mmty.N), M(mmty.M), p_data(mmty.p_data)
//    {}
    explicit Matrix (std::vector<T> va): name("Matrix "), N(va.size()), M(1), p_data(new T [N])
    {
        auto i_data = va.crbegin();
        for (T* p_curr = p_data + N - 1; p_curr >= p_data; --p_curr, ++i_data)
        {
            *p_curr = *i_data;
        }
        increase_name_count();
    }
    void edit_col (const int M_new)
    {
        if (!p_data)
        {
            p_data = new T [(N=1) * (M=M_new)];
        }
        T* p_curr = p_data;
        T* p_new = new T[N * M_new];
        T* p_end = p_new + N * (M_new - 1) - 1;
        --p_curr;
        for (T* p_curr_new = p_new - 1; p_curr_new <= p_end; p_curr += (M_new >= M) ? 0 : M - M_new)
        {
            for (int j = 0; j < M_new; ++j)
            {
                if (j < M)
                    *(++p_curr_new) = *(++p_curr);
                else
                {
                    for (; j < M_new; ++j)
                        *(++p_curr_new) = 0;
                #ifdef NDEBUG
            goto AfterNullN;
#endif // NDEBUG
                }
            }
        }
        AfterNullN:
        M = M_new;
        delete[] p_data;
        p_data = p_new;
    }
    void edit_row (const int a)
    {
        if (!p_data)
        {
            p_data = new T [(N=a) * (M=1)];
        }
        T* p_new = new T [a * M];
        T* p_end_new = p_new + a * M - 1;
        T* p_end = p_data + N * M - 2;
        T* p_curr = p_data - 1;
        for (T* p_curr_new = p_new; p_curr_new <= p_end_new; ++p_curr_new)
        {
            if (p_curr <= p_end)
                *p_curr_new = *(++p_curr);
            else
            {
                for (; p_curr_new <= p_end_new; ++p_curr_new)
                    *p_curr_new = 0;
#ifdef NDEBUG
                goto AfterNullR;
#endif // NDEBUG
            }
        }
        AfterNullR:
        N = a;
        delete[] p_data;
        p_data = p_new;
    }
    bool if_symmetric () const
    {
        if (N != M)
            return false;
        int i = 1;
        T* p_column_end = p_data;
        p_column_end += M - 1;
        for (T* p_row = p_data + N*M-1 - 1, *p_column = p_row - (M - 1); p_row >= p_data; p_row -= i, p_column += (N - i) * M - 1)
        {
            for (; p_column >= p_column_end; --p_row, p_column -= M)
            {
                if (*p_row != *p_column)
                {
                    return false;
                }
            }
            ++i;
        }
        return true;
    }
    T* get_pointer () const
    {
        return p_data;
    }
    int get_n () const
    {
        return N;
    }
    int get_m () const
    {
        return M;
    }
    int get_size () const
    {
        return N*M;
    }
    T* get_pa() const
    {
        return p_data;
    }
    Matrix operator+ (const Matrix& addend)
    {
        return plus_or_minus(addend, false);
    }
    Matrix operator- (const Matrix& subtrahend)
    {
        return plus_or_minus(subtrahend, true);
    }
    Matrix operator+ (const std::vector<T>& vect)
    {
        Matrix addend (vect);
        return *this+addend;
    }
    Matrix operator- (const std::vector<T>& vect)
    {
        Matrix subtrahend (vect);
        return *this-subtrahend;
    }
    Matrix& operator= (const Matrix& matrix)
    {
        if(p_data==matrix.p_data)
        {return *this;}
        delete[] p_data;
        name = matrix.name;
        N = matrix.N;
        M = matrix.M;
        p_data = new T [N * M];
        T* p_curr_new = matrix.p_data + N * M;
        for (T* p_curr = p_data + N * M - 1; p_curr >= p_data; --p_curr)
            *p_curr = *(--p_curr_new);
        return *this;
    }
    Matrix operator* (const Matrix& multiplier) const
    {
        assert(M == multiplier.N);
        T* p_new = new T [N * multiplier.M];
        for (T* ptr = p_new + N * multiplier.M - 1; ptr >= p_new; --ptr)
            *ptr = 0;
        T* p_end = p_data + N * M;
        T* p_curr_new = p_new;
        for (T* p_temp_curr = p_data; p_temp_curr < p_end; p_temp_curr += M)
            for (T* p_temp_m_curr = multiplier.p_data; p_temp_m_curr < multiplier.p_data + multiplier.M; ++p_temp_m_curr, ++p_curr_new)
                for (T* p_curr = p_temp_curr, *p_m_curr = p_temp_m_curr; p_curr < p_temp_curr + M; p_m_curr += multiplier.M, ++p_curr)
                    *p_curr_new += *p_curr * *p_m_curr;
        Matrix<T> answer (p_new, N, multiplier.M);
        delete[] p_new;
        return answer;
    }
    Matrix<T> multiple_matrix_by_number (const T num) const
    {
        Matrix<T> answer (N, M);
        for (T* p_curr_new = answer.p_data + answer.N * answer.M - 1, *p_curr = p_data + N * M - 1; p_curr_new >= answer.p_data; --p_curr_new, --p_curr)
        {
            *p_curr_new = *p_curr * num;
        }
        return answer;
    }
    explicit operator std::vector<T> ()
    {
        std::vector<T> va(0);
        if (M != 1)
        {
            std::cout << "Error while trynig to create vector from " << name << " with " << N << '*' << M << "dimensions";
            return va;
        }
        T* p_end = p_data + N;
        for (T* p_curr = p_data; p_curr < p_end; ++p_curr)
        {
            va.push_back(*p_curr);
        }
        return va;
    }
    T max_element() const
    {
        assert (p_data != nullptr);
        T* el = p_data + N * M - 1;
        for (T* p_curr = p_data + N * M - 2; p_curr >= p_data; --p_curr)
        {
            if (*el < *p_curr)
                *el = *p_curr;
        }
        return *el;
    }
    T max_element_abs() const
    {
        assert (p_data != nullptr);
        T el = fabs(*(p_data + N * M - 1));
        for (T* p_curr = p_data + N * M - 2; p_curr >= p_data; --p_curr)
        {
            if (el < fabs(*p_curr))
                el = fabs(*p_curr);
        }
        return el;
    }
    T min_element() const
    {
        assert (p_data != nullptr);
        T* el = p_data + N * M - 1;
        for (T* p_curr = p_data + N * M - 2; p_curr >= p_data; --p_curr)
        {
            if (*el > *p_curr)
                *el = *p_curr;
        }
        return *el;
    }
    T min_element_abs() const
    {
        assert (p_data != nullptr);
        T el = fabs(*(p_data + N * M - 1));
        for (T* p_curr = p_data + N * M - 2; p_curr >= p_data; --p_curr)
        {
            if (el > fabs(*p_curr))
                el = fabs(*p_curr);
        }
        return el;
    }
    virtual T& unsafe_index (const int i)
    {
        return *(p_data + i);
    }
    virtual T unsafe_index_c (const int i) const
    {
        return *(p_data + i);
    }
    virtual T& operator[] (const int i)
    {
        assert(i < N*M && i >= 0);
            return *(p_data + i);
    }
    virtual T operator[] (const int i) const
    {
        assert(i < N*M && i >= 0);
            return *(p_data + i);
    }
#ifndef NDEBUG
    class iterator
    {
        T* p;
        T* const p_data;
        int size;
    public:
        Matrix<T>::iterator& operator= (const Matrix<T>::iterator& a)//TODO Why? !!! PRIVATE !!!
        {
            if(this==&a)
                return *this;
            assert (p_data == a.p_data && size == a.size);
            p = a.p;
            return *this;
        }
    public:
        iterator (T* a, T* const pb, const int s): p(a), p_data(pb), size(s)
        {
            assert (a < pb + s);
        }
        iterator (): p(nullptr), p_data(nullptr), size(0)
        {}
        Matrix<T>::iterator& operator++ ()
        {
            assert (p <= p_data + size - 1);
#ifndef NDEBUG
            if (p == p_data + size - 1)
                std::cout<<"Warning! p == p_data  + size - 1 in prefix ++"<<std::endl;
#endif // NDEBUG
            ++p;
            return *this;
        }
        const Matrix<T>::iterator operator++ (int)
        {
            assert (p <= p_data + size - 1);
#ifndef NDEBUG
            if (p == p_data + size - 1)
                std::cout<<"Warning! p == p_data  + size - 1 in postfix ++"<<std::endl;
#endif // NDEBUG
            return Matrix<T>::iterator(p++, p_data, size);
        }
        Matrix<T>::iterator& operator-- ()
        {
            assert (p >= p_data);
#ifndef NDEBUG
            if (p == p_data)
                std::cout<<"Warning! p == p_data in prefix --"<<std::endl;
#endif // NDEBUG
            --p;
            return *this;
        }
        const Matrix<T>::iterator operator-- (int)
        {
            assert (p >= p_data);
#ifndef NDEBUG
            if (p == p_data)
                std::cout<<"Warning! p == p_data in postfix --"<<std::endl;
#endif // NDEBUG
            return Matrix<T>::iterator(p--, p_data, size);
        }
        Matrix<T>::iterator operator-= (int i)
        {
            assert (p >= p_data);
#ifndef NDEBUG
            if (p - i < p_data)
                std::cout<<"Warning! p < p_data in -="<<i<<std::endl;
#endif // NDEBUG
            p -= i;
            return *this;
        }
        Matrix<T>::iterator operator+= (int i)
        {
            assert (p < p_data + size);
#ifndef NDEBUG
            if (p + i >= p_data + size)
                std::cout<<"Warning! p + i >= p_data + size in +="<<i<<std::endl;
#endif // NDEBUG
            p += i;
            return *this;
        }
        Matrix<T>::iterator operator+ (const int i) const
        {
            assert (p < p_data + size);
#ifndef NDEBUG
            if (p + i >= p_data + size)
                std::cout<<"Warning! p + i >= p_data + size in + "<<i<<std::endl;
#endif // NDEBUG
            Matrix<T>::iterator a = *this;
            a.p += i;
            return a;
        }
        Matrix<T>::iterator operator- (const int i) const
        {
            assert (p >= p_data);
#ifndef NDEBUG
            if (p - i < p_data)
                std::cout<<"Warning! p < p_data in - "<<i<<std::endl;
#endif // NDEBUG
            Matrix<T>::iterator a = *this;
            a.p -= i;
            return a;
        }
        T& operator* ()
        {
            assert (p != nullptr);
            return *p;
        }
        T operator* () const
        {
            assert (p != nullptr);
            return *p;
        }
        operator T* ()
        {
            return p;
        }
        bool operator== (const Matrix<T>::iterator& a) const
        {
            assert (p_data == a.p_data);
            return this->p == a;
        }
        bool operator<= (const Matrix<T>::iterator& a) const
        {
            assert (p_data == a.p_data);
            return p <= a.p;
        }
        bool operator>= (const Matrix<T>::iterator& a) const
        {
            assert (p_data == a.p_data);
            return p >= a.p;
        }
        bool operator> (const Matrix<T>::iterator& a) const
        {
            assert (p_data == a.p_data);
            return p > a.p;
        }
        bool operator< (const Matrix<T>::iterator& a) const
        {
            assert (p_data == a.p_data);
            return p < a.p;
        }
    };
//#else
//    using Matrix<T>::iterator T*
#endif // NDEBUG
#ifndef NDEBUG
    Matrix::iterator first_i () const
    {
        return Matrix::iterator(p_data, p_data, N*M);
    }
    Matrix::iterator last_i () const
    {
        return Matrix::iterator(p_data + get_size() - 1, p_data, N*M);
    }
#else // NDEBUG
    T* first_i () const
    {
        return p_data;
    }
    T* last_i () const
    {
        return p_data + get_size() - 1;
    }
#endif // NDEBUG
    virtual const Matrix<T>& view (const int bytes_per_element = 8, const bool show_name = false) const
    {
        int n_str = N;
        T* p_curr = p_data;
        if (show_name)
            std::cout << name << std::endl;
        while (n_str--)
        {
            for (int j = 0; j < M - 1; j++)
            {
                std::cout << std::setw(bytes_per_element) << (*p_curr) << " ";
                p_curr++;
            }
            std::cout<<std::setw(bytes_per_element)<<*p_curr++;
            std::cout<<std::endl;
        }
        return *this;
    }
    virtual bool read_from_file (std::ifstream& sa)
    {
        if (sa.eof())
            return false;
        sa.read(reinterpret_cast<char*>(&N), sizeof(int));
        if (sa.eof())
            return false;
        sa.read(reinterpret_cast<char*>(&M), sizeof(int));
        p_data = new T [N * M];
        if (!p_data)
            return false;
        const T* p_end = p_data + N * M;
        const int size = sizeof(T);
        for (T* p_curr = p_data; p_curr < p_end; ++p_curr)
        {
            if (sa.eof())
                return false;
            sa.read(reinterpret_cast<char*>(p_curr), size);
        }
        return true;
    }
    virtual void write_to_file (std::ofstream& sa) const
    {
        //string a = name;
        //a.erase(0, 7);
        //int Num = a;
        //sa.write(reinterpret_cast<char*>(), sizeof(int));
        sa.write(reinterpret_cast<const char*>(&N), sizeof(int));
        sa.write(reinterpret_cast<const char*>(&M), sizeof(int));
        const int size = sizeof(T);
        const T* p_end = p_data + N * M;
        for (T* p_curr = p_data; p_curr < p_end; ++p_curr)
        {
            sa.write(reinterpret_cast<char*>(p_curr), size);
        }
    }
    void fill_nulls ()
    {
        for (T* p_curr = p_data + N * M - 1; p_curr >= p_data; --p_curr)
            *p_curr = 0;
    }
    double norm_by_infinity () const
    {
        double* p_end = p_data + N * M;
        double answer = *p_data;
        for (double* p_curr = p_data + 1; p_curr < p_end; ++p_curr)
        {
            if (fabs(*p_curr) > fabs(answer))
                answer = *p_curr;
        }
        return answer;
    }
    double norm_by_euclid () const
    {
        double answer = 0;
        for (double* p_curr = p_data + N * M - 1; p_curr >= p_data; --p_curr)
        {
            answer += *p_curr * *p_curr;
        }
        return sqrt(answer);
    }
    double norm_by_max () const
    {
        double answer = 0;
        for (double* p_curr = p_data + N * M - 1; p_curr >= p_data; --p_curr)
        {
            answer += fabs(*p_curr);
        }
        return answer;
    }
    virtual ~Matrix()
    {
        delete[] p_data;
    }
};
template <class T>
int Matrix<T>::count = 0;





template <>
bool Matrix<double>::if_symmetric () const// TODO
{
    if (N != M)
        return false;
    int i = 1;
    double* p_column_end = p_data;
    for (double* p_row = p_data + N * M - 1 - 1, *p_column = p_row - (M - 1); p_row >= p_data; p_row -= i, p_column += (N - i) * M - 1)
    {
        for (; p_column >= p_column_end; --p_row, p_column -= M)
        {
            if (*p_row - *p_column > DBL_EPSILON)
            {
                return false;
            }
        }
        ++i;
    }
    return true;
}

template <>
bool Matrix<float>::if_symmetric () const
{
    if (N != M)
        return false;
    int i = 1;
    float* p_column_end = p_data;
    for (float* p_row = p_data +N*M-1 - 1, *p_column = p_row - (M - 1); p_row >= p_data; p_row -= i, p_column += (N - i) * M - 1)
    {
        for (; p_column >= p_column_end; --p_row, p_column -= M)
        {
            if (*p_row - *p_column > DBL_EPSILON)
            {
                return false;
            }
        }
        ++i;
    }
    return true;
}


//class Array_of_Functions2
//{
//    int N;
//    Function2* Fn;
//public:
//    Array_of_Functions2 (void): N(1), Fn (new Function2[N])
//    {}
//    Array_of_Functions2 (const int S): N(S), Fn (new Function2[N])
//    {}
//    Array_of_Functions2 (Array_of_Functions2& Ar): N(Ar.N), Fn (new Function2[N])
//    {
//        for (Function2* pFn = Fn + N - 1, *pFa = Ar.Fn +  Ar.N - 1; pFn >= Fn; --pFn, --pFa)
//        {
//            *pFn = *pFa;
//        }
//    }
//    Function2* Get_Fn (void)
//    {
//        return Fn;
//    }
//    Function2 operator[] (const int S) const
//    {
//        if (S < N)
//            return NullF2;
//        return *(Fn + S);
//    }
//    Matrix<double> f (const Matrix<double> X) const
//    {
//        Matrix<double> A (N, 1);
//        double* p_data = A.get_pa() + N - 1;
//        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn, --p_data)
//        {
//            *p_data = pFn->f(X);
//        }
//        return A;
//    }
//    Matrix_SLE Derivative (const Matrix<double> X) const
//    {
//        int Num_of_variables = X.get_size();
//        Matrix<double> A (N, Num_of_variables);
//        double* p_data = A.get_pa() + A.get_size() - 1;
//        Matrix<double> Temp (Num_of_variables, 1);
//        double* pTemp = 0;
//        for (Function2* pFn = Fn + N*Num_of_variables - 1; pFn >= Fn; --pFn, --p_data)
//        {
//            Temp = pFn->Derivative_by_definition_M_in_column(X);
//            pTemp = Temp.get_pa() + Num_of_variables - 1;
//            for (int i = 0; i < Num_of_variables; ++i, --p_data, --pTemp)
//            {
//                *p_data = *pTemp;
//            }
//        }
//        return A;
//    }
//    Matrix<double> Integral (const double from_x, const double from_y, const double to_x, const double to_y) const;
//    double Max (const Matrix<double> X)
//    {
//        double m = 0;
//        double v = 0;
//        for (Function2* pA = Fn + N - 1; pA >= Fn; --pA)
//        {
//            if (m < (v = pA->f(X)))
//            m = v;
//        }
//            return m;
//    }
////    ~Array_of_Functions2 (void)
////    {
////        delete[] Fn;
////    }
//};











//using Fune = double (*)(double);
//using pFune = double (**) (double);
//class Matrix_of_Functions
//{
//    int N;
//    double (**Fn)(double);
//public:
//    Matrix_of_Functions (void): N(1), Fn (new Fune [N])
//    {}
//    Matrix_of_Functions (const int S): N(S), Fn (new Fune[N])
//    {}
//    pFune Get_Fn (void)
//    {
//        return Fn;
//    }
//    Fune operator[] (const int S)
//    {
//        if (S < N)
//        return *(Fn + S);
//        return NullF;
//    }
//    Matrix_SLE()
//    };

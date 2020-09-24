#ifndef MATH_H
#define MATH_H
#include "../Lib/Matrix.h"
#include "Different.h"

double norm (std::vector<double> vect)
{
    Matrix<double> A(std::move(vect));
    double* p_end = A.data() + A.get_n() * A.get_m();
    double answer = *A.data();
    for (double* p_curr = A.data() + 1; p_curr < p_end; ++p_curr)
    {
        if (fabs(*p_curr) > fabs(answer))
            answer = *p_curr;
    }
    return answer;
}

class Matrix_SE: public Matrix<double>
{
    bool high_accuracy;
protected:
    char type;
public:
    Matrix_SE  (): Matrix<double>(), high_accuracy(true), type('N')
    {}
    Matrix_SE (const Matrix<double>& matrix): Matrix<double>(matrix), high_accuracy(true), type('N')
    {}
    Matrix_SE (double* ptr, const size_t n, const size_t m): Matrix<double>(ptr, n, m), high_accuracy(true), type('N')
    {}
    Matrix_SE (const Matrix_SE& matrix): Matrix<double>(matrix), high_accuracy (matrix.high_accuracy), type('N')
    {}
    //    Matrix_SE(Matrix_SE&& mmty) noexcept: high_accuracy(mmty.high_accuracy), type(mmty.type)
//    {}
    Matrix_SE (const size_t n, const size_t m): Matrix<double> (n, m), high_accuracy(true), type('N')
    {}
    void set_high_accuracy (const bool a)
    {
        high_accuracy = a;   }
    bool get_high_accuracy () const
    {
        return high_accuracy;
    }
    void multi_row_by_num (const size_t number_of_row, const double num)
    {
        double* p_begin = p_data + number_of_row * M;
        double* p_end = p_begin + M;
        for (double* p_curr = p_begin; p_curr < p_end; ++p_curr)
        {
            *p_curr = num * (*p_curr);
        }
    }
    void one_minus_second_row (const size_t num_of_one, const size_t num_of_two, const double num)
    {
        double* p_1 = p_data + num_of_one * M;
        const double* const p_end_1 = p_1 + M;
        for (double* p_2 = p_data + num_of_two * M; p_1 < p_end_1; ++p_1, ++p_2)
        {
            *p_1 = *p_1 - num * *p_2;
        }
    }
    void replace_rows (const size_t num_of_one, const size_t num_of_two)
    {
        double* p_end_1 = p_data + (num_of_one + 1) * M;
        double* p_end_2 = p_data + (num_of_two + 1) * M;
        if (!high_accuracy)
            for (double *p_1 = p_end_1 - M, *p_2 = p_end_2 - M; p_1 < p_end_1; ++p_1, ++p_2)
            {
                *p_1 = *p_1 + *p_2;
                *p_2 = *p_1 - *p_2;
                *p_1 = *p_1 - *p_2;
            }
        else
        {
            double temp = 0;
            for (double *p_1 = p_end_1 - M, *p_2 = p_end_2 - M; p_1 < p_end_1; ++p_1, ++p_2)
            {
                temp = *p_1;
                *p_1 = *p_2;
                *p_2 = temp;
            }
        }
    }
    size_t find_max_not_null_el_in_col (const size_t number_of_col) const
    {
        double* p_curr = p_data + number_of_col * (M + 1);
        const double* const p_end = p_data + N * M;
        double el = *p_curr;
        size_t index_of_max = number_of_col;
        size_t curr_index = number_of_col;
        while ((p_curr += M) < p_end)
        {
            curr_index++;
            if (*p_curr > el && fabs(*p_curr) > DBL_EPSILON)
            {
                el = *p_curr;
                index_of_max = curr_index;
            }
        }
        return index_of_max;
    }
    const Matrix_SE& view () const override
    {
        Matrix::view(8, false);
        return *this;
    }
    virtual bool read_from_file (const char* str)
    {
        std::ifstream ifsa (str, std::ios::binary);
        bool temp = read_from_file(ifsa);
        ifsa.close();
        return temp;
    }
    bool read_from_file (std::ifstream& sa) override
    {
        if (sa.eof())
            return false;
        sa.read(&type, sizeof(char));
        switch (type)
        {
            case 'N':
                break;
            case 'L':
                break;
            default:
                std::cerr<<"Attempt to read Matrix_SLE from not-Matrix_SLE file";
                return false;
        }
        if (sa.eof())
            return false;
        sa.read(reinterpret_cast<char*>(&N), sizeof(N));
        if (sa.eof())
            return false;
        sa.read(reinterpret_cast<char*>(&M), sizeof(M));
        delete[] p_data;
        p_data = new double [N * M];
        const double* p_end = p_data + N * M;
        const int size = sizeof(double);
        for (double* p_curr = p_data; p_curr < p_end; ++p_curr)
        {
            if (sa.eof())
                return false;
            sa.read(reinterpret_cast<char*>(p_curr), size);
        }
        sa.close();
        return true;
    }
    virtual void write_to_file (const char* str) const
    {
        std::ofstream ofsa (str, std::ios::binary);
        write_to_file(ofsa);
        ofsa.close();
    }
    void write_to_file (std::ofstream& sa) const override
    {
        const double* p_end = p_data + N * M;
        sa.write (&type, sizeof(char));
        sa.write(reinterpret_cast<const char*>(&N), sizeof(N));
        sa.write(reinterpret_cast<const char*>(&M), sizeof(M));
        const int size = sizeof(double);
        for (double* p_curr = p_data; p_curr < p_end; ++p_curr)
        {
            sa.write(reinterpret_cast<char*>(p_curr), size);
        }
        sa.close();
    }
};

class Matrix_SLE: public Matrix_SE
{
public:
    Matrix_SLE (): Matrix_SE()
    {
        type = 'L';
    }
    Matrix_SLE (const Matrix<double>& matrix): Matrix_SE(matrix)
    {
        type = 'L';
    }
    Matrix_SLE (double* ptr, const size_t a, const size_t b): Matrix_SE(ptr, a, b)
    {
        type = 'L';
    }
    Matrix_SLE (const Matrix_SLE& matrix_sle): Matrix_SE(matrix_sle)
    {
        type = 'L';
    }
//    Matrix_SLE(Matrix_SLE&& A) noexcept{}
    Matrix_SLE (const size_t n, const size_t m): Matrix_SE(n, m)
    {
        type = 'L';
    }
    Matrix_SLE T () const
    {
        Matrix_SLE answer (M, N);
        double* p_ans = answer.p_data;
        double* p_end_ans = p_ans + N * M;
        double* p_end = p_data + N * M;
        for (double* p_row = p_data, *p_col_ans = p_ans; p_row < p_end; p_row+=M, ++p_col_ans)
        {
            for (double* p_curr = p_row, *p_curr_ans = p_col_ans; p_curr_ans < p_end_ans; ++p_curr, p_curr_ans+=N)
            {
                *p_curr_ans = *p_curr;
            }
        }
        return answer;
    }
    Matrix<double> solve (const Matrix<double>& b) const
    {
        return (if_symmetric()) ? LDLT(b) : solve_by_Gauss_method(b);
    }
    Matrix<double> solve_by_Gauss_method (Matrix_SE b) const
    {
        assert(N==M);
        assert(N==b.get_size());
        Matrix_SLE work_matrix (*(this));
        size_t index = 0;
        double koef = 0;
        for (size_t i = 0; i < N - 1; ++i)
        {
            index = work_matrix.find_max_not_null_el_in_col(i);
            if (fabs(work_matrix[index * M + index]) < 10 * DBL_EPSILON)
            {
                std::cerr<<"System";
                this->view();
                std::cerr<<"haven't solutions"<<std::endl;
                assert (fabs(work_matrix[index * M + index]) < 10 * DBL_EPSILON);
            }
            if (i != index)
            {
                work_matrix.replace_rows(i, index);
                b.replace_rows(i, index);
            }
            koef = work_matrix[i * M + i];
            work_matrix.multi_row_by_num(i, 1 / koef);
            b.multi_row_by_num(i, 1 / koef);
            for (index = i + 1; index < N; ++index)
            {
                koef = work_matrix[index * M + i];
                work_matrix.one_minus_second_row(index, i, koef);
                b.one_minus_second_row(index, i, koef);
            }
        }
        Matrix<double> answer (N, 1);
        const double* p_matrix = work_matrix.get_ptr() + work_matrix.get_n() * work_matrix.get_m() - 1;
        const double* p_b = b.get_ptr() + b.get_n() * 1 - 1;
        double* const p_beg = answer.get_ptr() + answer.get_n() * 1 - 1; //const_?
        const double* const p_end = answer.get_ptr();
        //assert(index == N - 1);
        index = M;
        for (double* p_curr_var = p_beg; p_curr_var >= p_end; --p_curr_var, p_matrix -= (index--))
        {
            *p_curr_var = *(p_b--);
            for (double* p_curr_other_var = p_beg; p_curr_other_var > p_curr_var; --p_matrix, --p_curr_other_var)
                *p_curr_var = *p_curr_var - *p_matrix * (*p_curr_other_var);
            *p_curr_var = *p_curr_var / (*p_matrix);
        }
        return answer;
    }
    Matrix<double> LDLT (const Matrix<double>& b) const//TODO
    {
        Matrix_SLE L (N, M);
        Matrix_SLE D (N, M);
        Matrix_SLE L_T (N, M);
        double* pl_end = L.get_ptr() + N * M;
        double* pl = L.get_ptr();
        double* pd = D.get_ptr();
        for (size_t i = 0; i < N - 1; )
        {
            pd[i*M + i] = p_data[i * M + i];
            for (size_t k_1 = 0; k_1 < i; ++k_1)
            {
                pd[i*M + i] -= L[i*M + k_1] * L[i * M + k_1] * pd[k_1 * M + k_1];
            }
            ++i;
            for (size_t k = 0; k < i; ++k)
            {
                L[i*M+k] = p_data[i * M + k];
                for (size_t t = 0; t < k; ++t)
                    L[i*M+k] -= L[i*M+t]*L[k*M+t]*pd[t*M+t];
                L[i*M+k] /= pd[k*M+k];
            }
        }
        size_t i = N - 1;
        pd[i*M + i] = p_data[i * M + i];
        for (size_t k_1 = 0; k_1 < i; ++k_1)
        {
            pd[i*M + i] -= L[i*M + k_1] * L[i * M + k_1] * pd[k_1 * M + k_1];
        }
        for (; pl < pl_end; pl+= M + 1)
        {
            *pl = 1;
        }
        L_T = L.T();
        Matrix<double> answer (L.get_n(), 1);
        answer[0] = b[0];
        auto p_l_row = L.begin() + M;
        auto p_l_end = L.begin() + M;
        for (auto pb = b.begin() + 1, p_ans_curr = answer.begin() + 1; pb <= b.end(); ++pb, ++p_ans_curr)
        {
            *p_ans_curr = *pb;
            for (auto p_l = p_l_row, p_ans_curr_1 = p_ans_curr - 1; p_l >= p_l_end; --p_l, --p_ans_curr_1)
            {
//                if (p_l == p_l_end)
//                    assert (p_ans_curr_1 == answer.begin());
//                assert (*p_l != 0);
//                assert (*p_l != 1);
                *p_ans_curr -= *p_l * *p_ans_curr_1;
            }
            p_l_end += get_m();
            p_l_row += get_m() + 1;
        }
        for (auto pD = D.begin(), pY = answer.begin(); pD <= D.end(); pD += M + 1, ++pY)
        {
//            assert (*pD != 0);
            *pY = *pY/(*pD);
        }
        auto p_LT_row = L_T.end() - M;
        auto p_LT_end = L_T.end() - M;
        for (auto pY = answer.end() - 1; pY >= answer.begin(); --pY)
        {
            for (auto pLT = p_LT_row, pY1 = pY + 1; pLT <= p_LT_end; ++pLT, ++pY1)
            {
//                if (pLT == p_LT_end)
//                    assert (pY1 == answer.end());
//                assert (*pLT != 0);
//                assert (*pLT != 1);
                *pY -= *pLT**pY1;
            }
            p_LT_end -= get_m();
            p_LT_row -= get_m() + 1;
        }
        return answer;
    }
    const Matrix_SLE& view () const override
    {
        Matrix::view(8, false);
        return *this;
    }
};

class Function
{
    double (*fn)(double);
    static double null_f (double)
    {
        return 0;
    }
public:
    explicit Function (double (*fn_1) (double)): fn (fn_1)
    {}
    Function (): fn (null_f)
    {}
    static double derivative_by_definition (double (*der) (double), const double x)
    {
        return der(x);
    }
    double derivative_by_definition (const double x) const
    {
        double eps = 0.1;
        return (fn(x + eps) - fn(x)) / eps;
    }
    double f(const double x) const
    {
        return fn(x);
    }
    static double integral_by_user_def (double (*anti_der) (double), const double from, const double to)
    {
        return anti_der(to) - anti_der(from);
    }
    double integral_by_Sympthon_method (const double from, const double to, const unsigned long long max_steps = 100000ULL) const
    {
        unsigned long long n_steps = 1ULL;
        double answer = 0;
        double length = to - from;
        while(length / n_steps > DBL_EPSILON && n_steps < max_steps)
        {
            n_steps *= 2ULL;
        }
        n_steps -=  2ULL;
        double h = length / n_steps;
        for (double x = from + h; x < to; x+=h)
        {
            answer += 4 * f(x);
            x+=h;
            answer += 2 * f(x);
        }
        answer += f(from);
        answer += f(to);
        answer *= h / 3;
        return answer;
    }
    double integral_by_Sympthon_method(const double from, const double to, const double epsilon) const
    {
        unsigned long long steps = 1000ULL;
        double val = 0;
        double prev_val = 0;
        double deflection = 15 * epsilon;
        while (deflection < fabs((val = integral_by_Sympthon_method(from, to, steps)) - prev_val))
        {
            steps = 2 * steps;
            prev_val = val;
        }
        return val;
    }
    double integral_by_trapeze_method (const double from, const double to, const unsigned long long max_steps = 100000ULL) const
    {
#ifndef LLONG_MAX
#define LLONG_MAX numeric_limits<long long>::max()
#endif
#ifndef ULLONG_MAX
#define ULLONG_MAX numeric_limits<unsigned long long>::max()
#endif
        unsigned long long num_of_steps = 1ULL;
        double integral = 0;
        double length = to - from;
        while(length / num_of_steps > DBL_EPSILON && num_of_steps < max_steps)
        {
            num_of_steps *= 2ULL;
        }
        num_of_steps -=  2ULL;
        double step = length / num_of_steps;
        for (double x = from + step; x < to; x += step, --num_of_steps)
        {
            integral += 2 * f(x);
        }
        if (num_of_steps > LLONG_MAX)
            num_of_steps = LLONG_MAX - num_of_steps;
        if (num_of_steps > 0 && num_of_steps < 100 * length)
        {
            while (--num_of_steps != ULLONG_MAX)
            {
                integral += 2 * f(to - std::numeric_limits<double>::epsilon());
            }
            ++num_of_steps;
        }
        assert (num_of_steps == 0);
        integral += f(from);
        integral += f(to);
        integral *= step / 2;
        return integral;
    }
    double integral_by_trapeze_method(const double from, const double to, const double epsilon) const
    {
        unsigned long long steps = 1000ULL;
        double val = 0;
        double prev_val = 0;
        double deflection = 3 * epsilon;
        while (deflection < fabs((val = integral_by_trapeze_method(from, to, steps)) - prev_val))
        {
            steps = 2 * steps;
            prev_val = val;
        }
        return val;
    }
};



/*class Matrix_SNE: public Matrix<Function>
{
    char Type;
public:
    Matrix_SNE (): Matrix<Function>(false), Type('N')
    {}
    Matrix_SNE (const size_t N1, const size_t M1): Matrix<Function>(N1, M1), Type('N')
    {}
    Matrix_SNE (Matrix_SNE& Ar): Matrix<Function> (*(reinterpret_cast< Matrix<Function>* > (&Ar))), Type ('N')
    {}
//    Function operator[] (const int S) const
//    {
//        if (S < N)
//            return null_f;
//        return *(fn + S);
//    }
    Matrix<double> f (const double x) const
    {
        Matrix<double> A (N, 1);
        double* pA = A.get_pa() + N - 1;
        for (Function* ppa = this->p_data + this->N - 1; ppa >= this->p_data; --(ppa), --pA)
        {
            *pA = ppa->f(x);
        }
        return A;
    }
//    Matrix_SLE derivative (const Matrix<double> X) const
//    {
//        int Num_of_variables = X.get_size();
//        Matrix<double> A (N, Num_of_variables);
//        double* pa = A.get_pa() + A.get_size() - 1;
//        Matrix<double> Temp (Num_of_variables, 1);
//        double* pTemp = 0;
//        for (Function_2* ppa = a + N*Num_of_variables - 1; ppa >= pa; --ppa, --pa)
//        {
//            Temp = ppa->derivative_by_definition_in_column(X);
//            pTemp = Temp.get_ptr() + Num_of_variables - 1;
//            for (int i = 0; i < Num_of_variables; ++i, --pa, --pTemp)
//            {
//                *pa = *pTemp;
//            }
//        }
//        return A;
//    }
//    Matrix<double> integral (const double from_x, const double from_y, const double to_x, const double to_y) const;
    double max_value (const double x)
    {
        double m = 0;
        double v = 0;
        for (Function* pA = p_data + N - 1; pA >= p_data; --pA)
        {
            if (m < (v = pA->f(x)))
                m = v;
        }
        return m;
    }
};*/
class Function_2
{
    double (*fn)(const Matrix<double>&);
    static double null_f (const Matrix<double>&)
    {
        return 0;
    }
public:
    explicit Function_2 (double (*fn_1) (const Matrix<double>& x)): fn (fn_1)
    {}
    Function_2 (): fn (null_f)
    {}
    static double derivative_by_definition (double (*der) (const Matrix<double>&), const Matrix<double>& x)
    {
        return der(x);
    }
    Matrix<double> derivative_by_definition_in_column (const Matrix<double>& X) const
    {
        size_t num_of_variables = X.get_size();
        Matrix<double> answer (num_of_variables, 1);
        double* p_ans = answer.get_ptr();
        const double M = 0.01; //0.05   0.1/////////
        double F_X = fn(X);
        Matrix<double> DX(X);
        double* p_end = DX.get_ptr() + DX.get_size();
        for (double* pDX = DX.get_ptr(), *pX = X.get_ptr(); pDX < p_end; ++pDX, ++pX, ++p_ans)
        {
            if( *pDX < 10*DBL_EPSILON)
            {
                *pDX += 1e-6;
                *p_ans = (fn(DX) - F_X) / 1e-6;
                *pDX = *pX;
            }
            else
            {
                *pDX += M*(*pDX);
                *p_ans = (fn(DX) - F_X) / (M * (*pX));
                *pDX = *pX;
            }
        }
        return answer;
    }
    double f(const Matrix<double>& X) const
    {
        return fn(X);
    }
    static double integral_by_definition (double (*anti_der) (const Matrix<double>&), const Matrix<double>& from, const Matrix<double>& to)
    {
        return anti_der(to) - anti_der(from);
    }
    /// By Sympthon method, only for 2 variables
    double integral_by_definition (const Matrix<double>& from, const Matrix<double>& to, const unsigned long long x_max_steps= 1000ULL, const unsigned long long y_max_steps= 1000ULL) const
    {
        assert (to.get_size() == 2 && from.get_size() == 2);
        double integral = 0;
        unsigned long long x_steps = 1;
        double x_length = to[0] - from[0];
        while(x_length / x_steps > DBL_EPSILON && x_steps < x_max_steps)
        {
            x_steps *= 2;
        }
        double x_step = x_length / x_steps;
        double x_begin = from[0];
        double x_end = to[0];
        unsigned long long y_steps = 1;
        double y_length = to[1] - from[1];
        while(y_length / y_steps > DBL_EPSILON && y_steps < y_max_steps)
        {
            y_steps *= 2;
        }
        double y_step = y_length / y_steps;
        double y_end = to[1];
        Matrix<double> temp;
        for (Matrix<double> curr_point = from; curr_point.at(1) < y_end; curr_point.at(1)+= y_step + y_step)
        {
            for (curr_point.at(0) = x_begin; curr_point.at(0) < x_end; curr_point.at(0)+= x_step + x_step)
            {
                temp = curr_point;
                integral += f(temp);
                temp.at(0) += x_step;
                integral += 4 * f(temp);
                temp.at(0) += x_step;
                integral += f(temp);

                temp.at(0) = curr_point.at(0);
                temp.at(1) += y_step;
                integral += 4 * f(temp);
                temp.at(0) += x_step;
                integral += 16 * f(temp);
                temp.at(0) += x_step;
                integral += 4 * f(temp);

                temp.at(0) = curr_point.at(0);
                temp.at(1) += y_step;
                integral += f(temp);
                temp.at(0) += x_step;
                integral += 4 * f(temp);
                temp.at(0) += x_step;
                integral += f(temp);
            }
        }
        return integral*= x_step * y_step / 9;
    }
};
class Array_of_functions_2
{
    size_t N;
    Function_2* fn;
public:
    Array_of_functions_2 (): N(1), fn (new Function_2[N])
    {}
    explicit Array_of_functions_2 (const size_t size): N(size), fn (new Function_2[N])
    {}
    Array_of_functions_2 (const Array_of_functions_2& arr): N(arr.N), fn (new Function_2[N])
    {
        for (Function_2* p_fn = fn + N - 1, *p_ar = arr.fn + arr.N - 1; p_fn >= fn; --p_fn, --p_ar)
        {
            *p_fn = *p_ar;
        }
    }
    Array_of_functions_2& operator= (const Array_of_functions_2& arr)
    {
        N = arr.N;
        if (fn == arr.fn)
        {
            return *this;
        }
        delete[] fn;
        fn = new Function_2[arr.get_size()];
        Function_2* p_arr = arr.fn + arr.N - 1;
        for (Function_2* p_curr = fn + N - 1; p_curr >= fn; --p_curr, --p_arr)
        {
            *p_curr = *p_arr;
        }
        return *this;
    }
    size_t get_size () const
    {
        return N;
    }
    Function_2* get_fn () const
    {
        return fn;
    }
    Function_2* data () const
    {
        return fn;
    }
    Function_2& operator[] (const size_t i) const
    {
        assert (i < N);
        return *(fn + i);
    }
    /// returns matrix in one column
    Matrix<double> f (const Matrix<double>& X) const
    {
        Matrix<double> values (N, 1);
        double* p_v = values.get_ptr() + N - 1;
        for (Function_2* p_curr = fn + N - 1; p_curr >= fn; --p_curr, --p_v)
        {
            *p_v = p_curr->f(X);
        }
        return values;
    }
    /// returns matrix in Num_of_variables column and N rows derivatives by definition
    Matrix<double> derivative (const Matrix<double>& X) const
    {
        const size_t num_of_variables = X.get_size();
        Matrix<double> answer (N, num_of_variables);
        double* p_curr = answer.get_ptr() + answer.get_size() - 1;
        Matrix<double> curr_der (num_of_variables, 1);
        double* p_curr_der = nullptr;
        for (Function_2* p_fn = fn + N - 1; p_fn >= fn; --p_fn)
        {
            curr_der = p_fn->derivative_by_definition_in_column(X);
            p_curr_der = curr_der.get_ptr() + num_of_variables - 1;
            for (size_t i = 0; i < num_of_variables; ++i, --p_curr, --p_curr_der)
            {
                *p_curr = *p_curr_der;
            }
        }
        return answer;
    }
    ///returns matrix in one column
    Matrix<double> integral (const Matrix<double>& from, const Matrix<double>& to) const
    {
        assert (from.get_size() == to.get_size() && (from.get_m() == to.get_m() || from.get_n() == to.get_n()));
        Matrix<double> answer (N, 1);
        double* p_curr = answer.end();
        for (Function_2* pf = fn + N - 1; pf == fn; --pf, --p_curr)
        {
            *p_curr = pf->integral_by_definition(from, to);
        }
        return answer;
    }
    double max_value (const Matrix<double>& X) const
    {
        double m = 0;
        double v = 0;
        for (Function_2* p_curr = fn + N - 1; p_curr >= fn; --p_curr)
        {
            if (m < (v = fabs(p_curr->f(X))))
                m = v;
        }
        return m;
    }
    ~Array_of_functions_2 ()
    {
        delete[] fn;
    }
};
extern std::string path;
namespace solve_nonlinear_equations
{
double f_delta_2 (const Matrix<double>& X_i, const Matrix<double>& delta)
{
    double delta_2 = 0;
    double diff = 0;
    for (double* p_x = X_i.get_ptr() + X_i.get_size() - 1, *dx = delta.get_ptr() + delta.get_size() - 1; p_x >= X_i.get_ptr(); --p_x, --dx)
    {
        if (fabs(*p_x + *dx) < 1 && fabs(*dx) > delta_2)
        {
            delta_2 = *dx;
        }
        else if ((diff = fabs(*dx / (*p_x + *dx))) > delta_2)
        {
            delta_2 = diff;
        }
    }
    return delta_2;
}
                    /// by Newton method
Matrix<double> solve_SNE (const Array_of_functions_2& func, const Array_of_functions_2& der_of_arr, const Matrix<double>& xy, const double eps_1, const double eps_2)
{
    assert (eps_1 > 1e-10);
    assert (eps_2 > 1e-10);
    const int num_it = 10000;
    int i = 0;
    Matrix<double> xy_i (xy);
    Matrix<double> delta(xy.get_n(), xy.get_m());
    Matrix<double> high_precision_xy_i (xy);
    Matrix<double> high_precision_delta(xy.get_n(), xy.get_m());
    Matrix_SLE high_precision_derivative (xy.get_n(), xy.get_n()); ///It is not a mistake!
    Matrix_SLE derivative;
    Matrix<double> solution;
    double delta_1 = 0;
    double delta_2 = 0;
    std::ofstream ofstream(path+"2_Lab.log");
    ofstream<<"With user defined"<<std::endl;
    do
    {
        xy_i = xy_i + delta;
        solution = func.f(xy_i).multiply_matrix_by_number(-1);
        delta_1 = solution[0];
        derivative = func.derivative(xy_i);
        delta = derivative.solve(solution);
        high_precision_xy_i = high_precision_xy_i + high_precision_delta;
        Function_2* pF = der_of_arr.get_fn() + der_of_arr.get_size() - 1;
        for (double* pD = high_precision_derivative.get_ptr() + high_precision_derivative.get_size() - 1; pD >=
                                                                                                          high_precision_derivative.get_ptr(); --pD, --pF)
        {
            *pD = -pF->f(high_precision_xy_i);// Minus against write solution = func.f(high_precision_xy_i)
        }                         //.multiply_matrix_by_number(-1);
        solution = func.f(high_precision_xy_i);
        high_precision_delta = high_precision_derivative.solve(solution);
        if (fabs(delta_1) > fabs(solution[0]))
        {
            ofstream<<"1 better"<<std::endl;
        }
        else
        {
            ofstream<<"0 worse"<<std::endl;
        }
        ++i;
    }
    while ((delta_1 = func.max_value(xy_i)) > eps_1 && (delta_2 = f_delta_2(xy_i, delta)) > eps_2 && i <= num_it);
    return xy_i;
}
}


namespace solve_differential_equations
{
    std::vector<Matrix<double> > explicit_Euler_method (const Array_of_functions_2& F, const double x_from, const double x_to, const Matrix<double>& u_0, const Matrix<double>& epsilon, const double max_step)
{
    assert (F.data() != nullptr);
    assert (u_0.data() != nullptr);
    assert (epsilon.data() != nullptr);
    assert (x_from < x_to);
    assert (F.get_size() == epsilon.get_size());
    assert (F.get_size() == u_0.get_size());
    assert (max_step > DBL_EPSILON);
    const size_t N = F.get_size();
    double x = x_from;
    const int iterations = 150000;
    Matrix<double> yk (N, iterations);
    Matrix<double> tk (1, iterations);
    Matrix<double> f_curr = u_0;
    Matrix<double> tk_curr (1, N);
    Matrix<double> y_curr (1, N+1);
    auto tk_i = tk.begin();
    *tk_i = x_from;
    ++tk_i;
    const auto eps_end = epsilon.end();
    auto yk_first = yk.begin();
    auto yk_end = yk.begin() + N;
    for (auto iu = u_0.begin() + (N - 1), iy = yk_end - 1; iu >= u_0.begin(); --iu, --iy)
    {
        *iy = *iu;
    }
    int i = 0;
    while (x < x_to && i < iterations)
    {
        for (auto iy = y_curr.begin(), iyk = yk_first; iyk < yk_end; ++iy, ++iyk)
        {
            *iy = *iyk;
        }
        *y_curr.end() = x;
        yk_first += N;
        yk_end += N;
        f_curr = F.f(y_curr);
        for (auto it = tk_curr.begin(), iF = f_curr.begin(), iEps = epsilon.begin(); iEps <= eps_end; ++it, ++iF, ++iEps)
        {
            *it = *iEps/(fabs(*iF)+ *iEps / max_step);
        }
        *tk_i = tk_curr.min_element();
        for (auto it = yk_first, iF = f_curr.begin(); it < yk_end; ++it, ++iF)
        {
            *it = *(it-N) + *tk_i**iF;
        }
        x += *tk_i;
        *tk_i = x;
        ++tk_i;
        ++i;
    }
        std::cout<<++i<<std::endl;
        std::vector<Matrix<double> > tk_and_yk (2);
    tk_and_yk[0] = tk;
    tk_and_yk[1]= yk;
    return tk_and_yk;
}
    std::vector<Matrix<double> > implicit_Euler_method (const Array_of_functions_2& F, const double x_from, const double x_to, const Matrix<double>& u_0, const Matrix<double>& epsilon, const double t_min, const double t_max, const bool if_use_three_zones_strategy = true)
{
    assert (F.data() != nullptr);
    assert (u_0.data() != nullptr);
    assert (epsilon.data() != nullptr);
    assert (F.get_size() == epsilon.get_size());
    assert (F.get_size() == u_0.get_size());
    assert (x_from < x_to);
    assert (t_min > DBL_EPSILON);
    assert (t_max > DBL_EPSILON);
    assert (t_min < t_max);
    for (auto i_e = epsilon.end(); i_e >= epsilon.begin(); --i_e)
        assert (*i_e > 0);
    const size_t N = F.get_size();
    int iterations = 20000;
    int it_curr = 0;
    double tk = t_min;
    double tk_prv = t_min;
    Matrix<double> tk_curr_next (1, N);
    Matrix<double> y_prv_curr (u_0); ///y on previous step
    Matrix<double> y_curr (u_0);
    Matrix<double> y_next_curr (u_0);///y on next step
    Matrix<double> x_curr (N + 1, 1);
    Matrix<double> epsilon_curr(epsilon.get_n(), epsilon.get_m());
    Matrix<double> yk (N, iterations);
    Matrix<double> tk_final (1, iterations);
    x_curr.at(N) = x_from;
    auto yk_first = yk.begin() + N;
    auto yk_last = yk_first;
    for (auto iyk = yk.begin(), iu_0 = u_0.begin(); iyk < yk_first; ++iyk, ++iu_0)
    {
        *iyk = *iu_0;
    }
    tk_final.at(0) = x_from;
    auto itk_final = tk_final.begin() + 1;
    const int num_of_it_in_Newton = 100;///number of iterations in Newton method
    int i = 0;
    const double eps_1 = 1e-4;
    const double eps_2 = 1e-4;
    double delta_1 = 0;
    double delta_2 = 0;
    Matrix<double> delta(u_0.get_n(), u_0.get_m());
    Matrix<double> f_curr;
    Matrix_SLE derivative;
    Matrix<double> arg_for_functions(2 * N + 2, 1);/// with x, t on k iteration
    auto i_arg = arg_for_functions.begin();
    while (x_curr.at(N) < x_to && it_curr < iterations)
    {

        x_curr.at(N) += tk;
        *itk_final = x_curr.at(N);
        ++itk_final;
        i_arg = arg_for_functions.end();
        *i_arg = tk;
        --i_arg;
        *i_arg = x_curr.at(N);
        --i_arg;
        for (auto iy_curr = y_curr.end(); iy_curr >= y_curr.begin(); --iy_curr, --i_arg)
        {
            *i_arg = *iy_curr;
        }
eps_curr_bigger_then_epsilon:
        do
        {
            y_next_curr = (y_next_curr + delta);
            i_arg = arg_for_functions.begin() + (N - 1);
            for (auto iy_next_curr = y_next_curr.end(); iy_next_curr >= y_next_curr.begin(); --iy_next_curr, --i_arg)
            {
                *i_arg = *iy_next_curr;
            }
            f_curr = F.f(arg_for_functions).multiply_matrix_by_number(-1.0);
            derivative = F.derivative(arg_for_functions);
            derivative.edit_col(N);
            delta = derivative.solve(f_curr);
        }
        while (++i <= num_of_it_in_Newton && ((delta_1 = f_curr.max_element_abs()) > eps_1) && ((delta_2 = solve_nonlinear_equations::f_delta_2(y_next_curr, delta)) > eps_2));
        delta.fill_nulls();
        i = 0;
        for (auto ie = epsilon_curr.end(), iy_prv = y_prv_curr.end(), iy_curr = y_curr.end(), iy_next = y_next_curr.end(); ie >= epsilon_curr.begin(); --ie, --iy_prv, --iy_curr, --iy_next)
        {
            *ie = tk/(tk+tk_prv)*((tk/tk_prv)*(*iy_curr - *iy_prv) + *iy_curr - *iy_next);
        }
        if (if_use_three_zones_strategy)
        {
            for (auto ie_curr = epsilon_curr.end(), ie = epsilon.end(); ie >= epsilon.begin(); --ie_curr, --ie)
            {
                if (*ie_curr > *ie)
                {
                    tk /= 2;
                    y_next_curr = y_curr;
                    x_curr.at(N) -= tk;
                    goto eps_curr_bigger_then_epsilon;
                }
            }
            for (auto itk_curr = tk_curr_next.end(), iEps_curr = epsilon_curr.end(), iEps = epsilon.end(); iEps >= epsilon.begin(); --itk_curr, --iEps, --iEps_curr)
            {
                *itk_curr = sqrt(*iEps/fabs(*iEps_curr))*tk;
            }
        }
        else
            for (auto itk_curr = tk_curr_next.end(), iEps_curr = epsilon_curr.end(), iEps = epsilon.end(); iEps >= epsilon.begin(); --itk_curr, --iEps, --iEps_curr)
            {
                *itk_curr = fabs(*iEps_curr) > *iEps ? tk/2 : fabs(*iEps_curr) > *iEps/4 ? tk : tk+tk;
            }
        tk_prv = tk;
        tk = tk_curr_next.min_element();
        if (tk > t_max)
            tk = t_max;
        y_prv_curr = y_curr;
        y_curr = y_next_curr;
        yk_last += N;
        for (auto iyk = yk_first, iy_curr = y_curr.begin(); iyk < yk_last; ++iyk, ++iy_curr)
        {
            *iyk = *iy_curr;
        }
        yk_first += N;
        ++it_curr;
    }
    std::cout << it_curr << std::endl;
    std::vector<Matrix<double> > tk_and_yk (2);
    tk_and_yk.at(0) = tk_final;
    tk_and_yk.at(1) = yk;
    return tk_and_yk;
}

}

namespace approximation
{
double dispersion (const Matrix<double>& X, const Matrix<double>& Y, const Matrix<double>& polinom, const unsigned int i)
{
    assert(polinom.get_size() == i + 1);
    assert(X.get_size() >= i + 1);
    double diff = 0;
    double answer = 0;
    for (auto pX = X.end(), pY = Y.end(); pY >= Y.data(); --pX, --pY)
    {
        diff = *pY;
        auto k = i;
        for (auto p_pol = polinom.end(); p_pol >= polinom.data(); --p_pol)
        {
            if (p_pol == polinom.data())
                assert (k == 0);
            diff -= *p_pol * Pow(*pX, k--);
        }
        answer += diff * diff;
    }
    answer=answer/(X.get_size() - i - 1);
    /// answer=sqrt(answer);
    return answer;
}
/// @returns p[0]+p[1]*x+p[2]*x^2...+p[n-1]*x^{n-1}
Matrix<double> find_polinom_m_power (const Matrix<double>& X, const Matrix<double>& Y, double temp)
{
    assert (X.get_size() == Y.get_size() && ((X.get_n() == 1 && Y.get_n() == 1) || (X.get_m() == 1 && Y.get_m() == 1)));
    if (temp < DBL_EPSILON || temp > X.get_size())
        temp = X.get_size();
    const unsigned int M = temp;
    unsigned int i = 2*M;
    double* sums = new double [2 * M + 1];
    sums[0] = 1;
    for (double* p_sum = sums + 2 * M; p_sum > sums; --p_sum)
    {
        temp = 0;
        for(const double* p_x = X.get_ptr() + X.get_size() - 1; p_x >= X.get_ptr(); --p_x)
        {
            temp += Pow(*p_x, i);
        }
        --i;
        *p_sum = temp;
    }
    Matrix_SLE X_pows (M + 1, M + 1);
    for (double* p_row = X_pows.end(); p_row >= X_pows.get_ptr(); p_row -= X_pows.get_m())
    {
        double* p_sum = sums + 2 * M - i++;
        for (double* p_curr = p_row; p_curr > p_row - M - 1; --p_curr, --p_sum)
            *p_curr = *p_sum;
    }
    delete[] sums;
    i = M;
    Matrix<double> y_sums (M + 1, 1);
    for (double* y_sum_curr = y_sums.get_ptr() + M; y_sum_curr > y_sums.get_ptr(); --y_sum_curr)
    {
        temp = 0;
        for(const double* p_arr = X.get_ptr() + X.get_size() - 1, *py = Y.get_ptr() + Y.get_size() - 1; p_arr >= X.get_ptr(); --p_arr, --py)
        {
            temp += Pow(*p_arr, i) * *py;
        }
        --i;
        *y_sum_curr = temp;
    }
    for (const double *py = Y.get_ptr() + Y.get_size() - 1; py >= Y.get_ptr(); --py)
    {
        *y_sums.get_ptr() += *py;
    }
    Matrix<double> parameters = X_pows.solve(y_sums);
    return parameters;
}
Matrix<double> find_polinom (const Matrix<double>& X, const Matrix<double>& Y)
{
    assert (X.get_size() == Y.get_size() && ((X.get_n() == 1 && Y.get_n() == 1) || (X.get_m() == 1 && Y.get_m() == 1)));
    Matrix<double> polinom;
    double* dispersions = new double [X.get_size()-1];
    double* p_curr = dispersions;
    for (int i = 1; i < X.get_size(); ++i)
    {
        polinom = find_polinom_m_power(X, Y, i);
        *(p_curr++) = dispersion(X, Y, polinom, i);
    }
    View(dispersions, X.get_size()-1, 1);
    size_t min_pow = X.get_size() - 1;///(X.get_size()-1) size, (-1) as (size-1), (+1) as it polinom koef.
    size_t curr_pow= min_pow - 1;///(-1) as next element already in min_pow.
    for (p_curr=dispersions+X.get_size()-3; p_curr >= dispersions; --p_curr, --curr_pow)
    {
        if(*p_curr < dispersions[min_pow - 1])
        {
            min_pow=curr_pow;
        }
    }
    delete[] dispersions;
    return find_polinom_m_power(X, Y, min_pow);
}
}
#endif // MATH_H

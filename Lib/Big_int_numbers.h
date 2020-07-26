//
// Created by alexander on 17/07/2020.
//

#ifndef NUM_METHODS_BIG_INT_NUMBERS
#define NUM_METHODS_BIG_INT_NUMBERS

#include "Different.h"
#include <boost/logic/tribool.hpp>
#include <cassert>
#include <vector>
#include <type_traits>
#include <iostream>

class Big_number
{
public:
    Big_number() : sign_plus(true), base(10)
    {
        number = std::vector<unsigned char>({0});
    }

    template<typename T>
    explicit Big_number(T num): sign_plus(num >= /*static_cast<T>*/(0)), number(std::vector<unsigned char>()), base(10)
    {
        static_assert(std::is_integral<T>::value, "You trying to construct number with not integral type");
        while (num)
        {
            number.push_back(num % base);
            num /= base;
        }
    }

    Big_number(const Big_number &big_number) : sign_plus(big_number.sign_plus), number(big_number.number), base(big_number.base)
    {}

    Big_number(Big_number &&big_number) noexcept: sign_plus(big_number.sign_plus), number(std::move(big_number.number)), base(big_number.base)
    {
        std::cout<<"lvl"<<std::is_lvalue_reference_v<Big_number><<std::endl;
        std::cout<<std::is_rvalue_reference_v<Big_number><<std::endl;
    }

    template<typename T>
    T Get_Number() const
    {
        static_assert(std::is_integral<T>::value, "You trying to construct number with not integral type");
        assert(sign_plus || std::is_signed<T>::value);
        assert(2ULL<<sizeof(T) >= Pow(base, number.size()));
//        static_assert(static_cast<T>(2) << (8 * sizeof(T)) >= number.size());
        T num = *number.crbegin();
#ifndef NDEBUG
        T num_for_safety = num;
#endif
        for (auto it = number.crbegin() + 1; it < number.crend(); ++it)
        {
            num = *it + 10 * num;
#ifndef NDEBUG
            assert(num >= num_for_safety);
            num_for_safety = num;
#endif
        }
        if (!sign_plus)
            num = -num;
        return num;
    }


    Big_number& operator=(Big_number big_number)
    {
        std::swap(*this, big_number);
        return *this;
    }

    Big_number operator+(const Big_number &big_number) const
    {
        return *this;
    }

    Big_number operator-(const Big_number &big_number) const
    {
        return *this;
    }

    Big_number operator*(const Big_number &big_number) const
    {
        return *this;
    }

    Big_number operator%(Big_number big_number) const
    {
        return Divide (big_number, MOD);
    }

    Big_number operator/(const Big_number &big_number) const
    {
        return Divide (big_number, DIVIDE);
    }

    Big_number operator+=(const Big_number &big_number) const
    {
        return *this;
    }

    Big_number operator-=(const Big_number &big_number) const
    {
        return *this;
    }

    Big_number operator*=(const Big_number &big_number) const
    {
        return *this;
    }

    Big_number operator/=(const Big_number &big_number) const
    {
        return *this;
    }

    boost::tribool operator==(const Big_number &big_number) const
    {
        return number == big_number.number && sign_plus == big_number.sign_plus && base == big_number.base;
    }

    boost::tribool operator!=(const Big_number &big_number) const
    {
        return !operator==(big_number);
    }

    boost::tribool operator<(const Big_number &big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    boost::tribool operator<=(const Big_number &big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    boost::tribool operator>(const Big_number &big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    boost::tribool operator>=(const Big_number &big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    ~Big_number()
    {

    }
private:
    const static bool MOD = false;
    const static bool DIVIDE = true;
    std::vector<unsigned char> number;
    ///@note Zero has +
    bool sign_plus;
    unsigned char base;
    Big_number Divide(const Big_number& divisor, const bool mode) const
    {
        return *this;
    }
};

#endif //NUM_METHODS_BIG_INT_NUMBERS

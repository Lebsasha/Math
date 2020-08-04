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
#include <algorithm>
#include <iostream>

class Big_number
{
public:
    Big_number() : /*sign_plus(true),*/ base(10)
    {
        number = std::vector<unsigned char>({0});
    }

    template<typename T>
    explicit Big_number(T num): /*sign_plus(num >= *//*static_cast<T>*//*(0)),*/ number(std::vector<unsigned char>()), base(10)
    {
        static_assert(std::is_integral<T>::value, "You trying to construct number with not integral type");
        while (num)
        {
            number.push_back(num % base);
            num /= base;
        }
        std::reverse(number.begin(), number.end());
    }

    Big_number(const Big_number& big_number) : /*sign_plus(big_number.sign_plus),*/ number(big_number.number), base(big_number.base)
    {
        std::cout << "I'm copy operator" << std::endl;
    }

    Big_number(Big_number&& big_number) noexcept: /*sign_plus(big_number.sign_plus),*/ number(std::move(big_number.number)),
                                                                                       base(big_number.base)
    {
        std::cout << "I'm move operator" << std::endl;
    }

    template<typename T>
    T Get_Number() const
    {
        static_assert(std::is_integral<T>::value, "You trying to construct number with not integral type");
//        assert(sign_plus || std::is_signed<T>::value);
        assert((2ULL << (sizeof(T) * 8 - 1)) - 1 >= Pow(base, number.size()) - 1);
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
//        if (!sign_plus)
//            num = -num;
        return num;
    }

    Big_number& operator=(Big_number big_number)// TODO swap error prone?
    {
        std::swap(*this, big_number);
        return *this;
    }

    void set_base(const unsigned char base_)
    {
        assert(base_ <= std::numeric_limits<unsigned char>::max() / 2);
        base = base_;

    }

    [[nodiscard]] unsigned char get_base() const
    {
        return base;
    }

    Big_number operator+(const Big_number& big_number) const
    {
        Big_number ans = this->number.size() >= big_number.number.size() ? *this : big_number;
        const Big_number& addend = this->number.size() >= big_number.number.size() ? big_number : *this;
        auto itEnd = addend.number.cend();
        auto digit_l = addend.number.cbegin();
        for (auto digit_r = ans.number.begin(); digit_l < itEnd; ++digit_r, ++digit_l)
        {
            *digit_r += *digit_l;
        }
        itEnd = ans.number.cend() - 1;
        for (auto digit_r = ans.number.begin(); digit_r < itEnd; ++digit_r)
        {
            *(digit_r + 1) += *digit_r / base;
            *digit_r = *digit_r % base;
        }
        if (*ans.number.rbegin() >= base)
        {
            ans.number.push_back(*ans.number.rbegin() / base);
            *(ans.number.rbegin() + 1) = *(ans.number.rbegin() + 1) % base;
        }

        return ans;
    }

    Big_number operator-(const Big_number& big_number) const
    {
        return *this;
    }

    Big_number operator*(const Big_number& big_number) const
    {
        Big_number ans;
        for(auto digit_l = big_number.number.crbegin(); digit_l <big_number.number.crend(); ++digit_l)
        {
//            ans += multiply_by_num (*digit_l);

        }
        return *this;
    }

    Big_number operator%(Big_number big_number) const
    {
        return Divide(big_number, MOD);
    }

    Big_number operator/(const Big_number& big_number) const
    {
        return Divide(big_number, DIVIDE);
    }

    Big_number& operator+=(const Big_number& big_number)
    {
        *this = *this + big_number;
        return *this;
    }

    Big_number& operator-=(const Big_number& big_number)
    {
        *this = *this - big_number;
        return *this;
    }

    Big_number& operator*=(const Big_number& big_number)
    {
        *this = *this * big_number;
        return *this;
    }

    Big_number& operator/=(const Big_number& big_number)
    {
        *this = *this / big_number;
        return *this;
    }

    bool operator==(const Big_number& big_number) const
    {
        return number == big_number.number/* && sign_plus == big_number.sign_plus*/ && base == big_number.base;
    }

    bool operator!=(const Big_number& big_number) const
    {
        return !operator==(big_number);
    }

    boost::tribool operator<(const Big_number& big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    boost::tribool operator<=(const Big_number& big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    boost::tribool operator>(const Big_number& big_number) const
    {
        return boost::tribool(boost::logic::indeterminate);
    }

    boost::tribool operator>=(const Big_number& big_number) const
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
    ///@note Zero has '+'
//    bool sign_plus;
    /// base can be in range [0, 127] on "standard" computer
    unsigned char base;

    Big_number Divide(const Big_number& divisor, const bool mode) const
    {
        return *this;
    }
};

#endif //NUM_METHODS_BIG_INT_NUMBERS

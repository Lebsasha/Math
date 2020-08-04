//
// Created by alexander on 17/07/2020.
//

#ifndef NUM_METHODS_BIG_INT_NUMBERS
#define NUM_METHODS_BIG_INT_NUMBERS

#include "Different.h"
#include <iterator>
//#include <cassert>
//#include <vector>
//#include <type_traits>
#include <algorithm>
//#include <iostream>

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
        T num = *number.crbegin();
        for (auto it = number.crbegin() + 1; it < number.crend(); ++it)
        {
            num = *it + 10 * num;
        }
//        if (!sign_plus)
//            num = -num;
        return num;
    }

    Big_number& operator=(Big_number big_number)
    {
        number = std::move(big_number.number);
        base = big_number.base;
//        std::swap(*this, big_number);
        return *this;
    }

    /// base can be in range [0, 16] on "standard" computer
    void set_base(const int base_)
    {
        assert(base_ <= sqrt(std::numeric_limits<unsigned char>::max() + 1));
        base = base_;
    }

    [[nodiscard]] int get_base() const
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
        ans.normalise();
        return ans;
    }

    Big_number operator-(const Big_number& big_number) const //TODO
    {
        return *this;
    }

    Big_number operator*(const Big_number& big_number) const
    {
        Big_number ans = multiply_by_num(*big_number.number.crbegin());
        auto itrEnd = ans.number.rend();
        for (auto digit_l = big_number.number.crbegin() + 1; digit_l < big_number.number.crend(); ++digit_l)
        {
            ans.number.push_back(*ans.number.rbegin());
            itrEnd = ans.number.rend() - 1;
            for (auto itr = ans.number.rbegin() + 1; itr < itrEnd; ++itr)
            {
                *itr = *(itr + 1);
            }
            *itrEnd = 0;
            ans += multiply_by_num(*digit_l);
        }
        return ans;
    }

    Big_number operator%(Big_number big_number) const //TODO
    {
        return Divide(big_number, MOD);
    }

    Big_number operator/(const Big_number& big_number) const //TODO
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

    bool operator<(const Big_number& big_number) const
    {
        if (base != big_number.base)
        {
            std::cerr << "base != big_number.base in operator<" << std::endl;
            return false;
        }
        if (number.size() != big_number.number.size())
            return number.size() < big_number.number.size();
        auto itr1 = number.crbegin();
        auto itr2 = big_number.number.crbegin();
        for (; itr2 < big_number.number.crend(); ++itr1, ++itr2)
        {
            if (*itr1 == *itr2)
                continue;
            return *itr1 < *itr2;
        }
        return false; // Numbers are equal
    }

    bool operator<=(const Big_number& big_number) const
    {
        return *this == big_number || *this < big_number;
    }

    bool operator>(const Big_number& big_number) const
    {
        if (base != big_number.base)
        {
            std::cerr << "base != big_number.base in operator>" << std::endl;
            return false;
        }
        if (number.size() != big_number.number.size())
            return number.size() > big_number.number.size();
        auto itr1 = number.crbegin();
        auto itr2 = big_number.number.crbegin();
        for (; itr2 < big_number.number.crend(); ++itr1, ++itr2)
        {
            if (*itr1 == *itr2)
                continue;
            return *itr1 > *itr2;
        }
        return false; // Numbers are equal
    }

    bool operator>=(const Big_number& big_number) const
    {
        return *this == big_number || *this > big_number;
    }

    void View() const
    {
        Big_number Temp = *this;
        Temp.set_base(10);
        std::reverse_copy(Temp.number.begin(), Temp.number.end(), std::ostream_iterator<int>(std::cout, ""));
        std::cout << std::endl;
    }

    ~Big_number() = default;

private:
    const static bool MOD = false;
    const static bool DIVIDE = true;
    std::vector<unsigned char> number;
//    ///@note Zero has '+'
//    bool sign_plus;
    /// base can be in range [0, 16] on "standard" computer
    unsigned char base;

    void normalise()
    {
        Big_number& ans = *this;
        auto itEnd = ans.number.cend() - 1;
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
    }

    [[nodiscard]] Big_number multiply_by_num(const unsigned char& i) const
    {
        assert(i < 16);
        Big_number ans = *this;
        for (auto& digit: ans.number)
            digit = digit * i;
        ans.normalise();
        return ans;
    }

    [[nodiscard]] Big_number Divide(const Big_number& divisor, const bool mode) const
    {
        return *this;
    }

};

#endif //NUM_METHODS_BIG_INT_NUMBERS

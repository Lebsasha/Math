//
// Created by alexander on 17/07/2020.
//

#ifndef BIG_INT_NUMBERS
#define BIG_INT_NUMBERS

#include <iterator>
#include <algorithm>
#include <ostream>
#include "Different.h"
// TODO Move to .cpp
class Big_number
{
public:
    Big_number() :number(std::vector<unsigned char>({0})), base(10)
    {}

    template<typename T>
    explicit Big_number(T num):number(std::vector<unsigned char>()), base(10)
    {
        static_assert(std::is_integral<T>::value, "You trying to construct number with not integral type");
        if (!num)
            number.push_back(num);
        while (num)
        {
            number.push_back(num % base);
            num /= base;
        }
    }

    Big_number(const Big_number& big_number) :number(big_number.number), base(big_number.base)
    {
#ifndef NDEBUG
        std::cout << "I'm copy operator" << std::endl;
#endif
    }

    Big_number(Big_number&& big_number) noexcept:number(std::move(big_number.number)), base(big_number.base)
    {
#ifndef NDEBUG
        std::cout << "I'm move operator" << std::endl;
#endif
    }

    Big_number& operator=(Big_number big_number)
    {
        number = std::move(big_number.number);
        base = big_number.base;
#ifndef NDEBUG
        std::cout << "I'm operator=" << std::endl;
#endif
        return *this;
    }

    template<typename T>
    [[nodiscard]] T get_number() const
    {
        static_assert(std::is_integral<T>::value, "You trying to construct number with not integral type");
        T num = 2;
        assert((num << (sizeof(T) * 8 - 1-std::is_signed_v<T>)) - 1 >= Pow(base, number.size()) - 1);
        num = *number.crbegin();
        for (auto it = number.crbegin() + 1; it < number.crend(); ++it)
        {
            num = *it + base * num;
        }
        return num;
    }

    /// base can be in range [2, 16] on "standard" computer
    void set_base(const unsigned char new_base)
    {
        assert(new_base >= 2);
        assert(new_base <= sqrt(std::numeric_limits<unsigned char>::max() + 1));
        Big_number ans;
        if (base != 10)
        {
            int i = 0;
            for (auto digit: number)
            {
                ans += Big_number(base).pow(i) * digit;
                ++i;
            }
            *this = ans;
        }
        if (new_base != 10)
        {
            ans.number.clear();
            while (static_cast<bool>(*this))
            {
                ans.number.push_back((*this % new_base).get_number<unsigned char>());
                *this /= Big_number(new_base);
            }
            *this = ans;
            base = new_base;
        }
    }

    [[nodiscard]] int get_base() const
    {
        return base;
    }

    Big_number operator+(const Big_number& big_number) const
    {
        /// addend_1 + addend_2 = sum;
        /// addend_1 == sum == ans here for optimisation purposes
        Big_number ans = this->number.size() >= big_number.number.size() ? *this : big_number;
        const Big_number& addend_2 = this->number.size() >= big_number.number.size() ? big_number : *this;
        auto it_end = addend_2.number.cend();
        auto digit_l = addend_2.number.cbegin();
        for (auto digit_r = ans.number.begin(); digit_l < it_end; ++digit_r, ++digit_l)
        {
            *digit_r += *digit_l;
        }
        ans.normalise();
        return ans;
    }

    template<typename T>
    Big_number operator+(const T& num) const
    {
        return *this + Big_number(num);
    }

    Big_number operator-(const Big_number& big_number) const
    {
        /// minuend - subtrahend = difference;
        /// minuend == difference == ans here for optimisation purposes
        assert (*this >= big_number);
        Big_number ans = *this;
        auto i_minuend = ans.number.begin();
        auto i_minuend_finder_non_zeros = i_minuend;
        for (auto v_subtrahend: big_number.number)
        {
            if (*i_minuend >= v_subtrahend)
                *i_minuend -= v_subtrahend;
            else
            {
                i_minuend_finder_non_zeros = i_minuend;
                while (*(++i_minuend_finder_non_zeros) == 0);
                *(i_minuend_finder_non_zeros) -= 1;
                for (auto p = i_minuend + 1; p < i_minuend_finder_non_zeros; ++p)
                    *p = base - 1;
                *i_minuend += base;
                *i_minuend -= v_subtrahend;
            }
            ++i_minuend;
        }
        for (auto i_lead_zeros = ans.number.rbegin();
             i_lead_zeros < ans.number.rend() - 1 && *i_lead_zeros == 0; i_lead_zeros = ans.number.rbegin())
        {
            ans.number.erase((i_lead_zeros + 1).base());
        }
        return ans;
    }

    template<typename T>
    Big_number operator-(const T& num) const
    {
        return *this - Big_number(num);
    }

    Big_number operator*(const Big_number& multiplier) const
    {
        /// multiplicand * multiplier = product;
        /// multiplicand == product == ans here for optimisation purposes
        if (*this == Big_number() || multiplier == Big_number())
            return Big_number();
        Big_number ans = multiply_by_num(*multiplier.number.crbegin());
        auto i_end = ans.number.rend();
        for (auto digit_l = multiplier.number.crbegin() + 1; digit_l < multiplier.number.crend(); ++digit_l)
        {
            ans.number.push_back(*ans.number.rbegin());
            i_end = ans.number.rend() - 1;
            for (auto i_mover = ans.number.rbegin() + 1; i_mover < i_end; ++i_mover)
            {
                *i_mover = *(i_mover + 1);
            }
            *i_end = 0;
            ans += multiply_by_num(*digit_l);
        }
        return ans;
    }

    template<typename T>
    Big_number operator*(const T& num) const
    {
        return *this * Big_number(num);
    }

    Big_number operator%(const Big_number& divisor) const
    {
        return divide(divisor, false);
    }

    template<typename T>
    Big_number operator%(const T& num) const
    {
        return *this % Big_number(num);
    }

    Big_number operator/(const Big_number& divisor) const
    {
        return divide(divisor, true);
    }

    template<typename T>
    Big_number operator/(const T& num) const
    {
        return *this / Big_number(num);
    }

    [[nodiscard]] Big_number pow(int i) const
    {
        assert(i >= 0);
        if (!i)
            return Big_number(1);
        Big_number ans = *this;
        while (--i)
            ans *= *this;
        return ans;
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
        return number == big_number.number&& base == big_number.base;
    }

    bool operator!=(const Big_number& big_number) const
    {
        return !operator==(big_number);
    }

    bool operator<(const Big_number& big_number) const
    {
        if (base != big_number.base && base != 10 && big_number.base != 10)
        {
            std::cerr << "base != big_number.base in operator< or operator>" << std::endl;
            return false;
        }
        if (number.size() != big_number.number.size())
            return number.size() < big_number.number.size();
        auto i_left_digit = number.crbegin();
        auto i_right_digit = big_number.number.crbegin();
        for (; i_right_digit < big_number.number.crend(); ++i_left_digit, ++i_right_digit)
        {
            if (*i_left_digit != *i_right_digit)
                return *i_left_digit < *i_right_digit;
        }
        return false; // Numbers are equal
    }

    bool operator<=(const Big_number& big_number) const
    {
        return *this == big_number || *this < big_number;
    }

    bool operator>(const Big_number& big_number) const
    {
        return big_number < *this;
    }

    bool operator>=(const Big_number& big_number) const
    {
        return *this == big_number || *this > big_number;
    }

    explicit operator bool() const
    {
        return !(number.size() == 1 && number[0] == 0);
    }

    void view(std::ostream& ostr = std::cout) const
    {
        if (base != 10)
        {
            Big_number temp = *this;
            temp.set_base(10);
            std::reverse_copy(temp.number.begin(), temp.number.end(), std::ostream_iterator<int>(ostr, ""));
        }
        else
        {
            const Big_number& temp = *this;
            std::reverse_copy(temp.number.begin(), temp.number.end(), std::ostream_iterator<int>(ostr, ""));
        }
        ostr << std::endl;
    }

    ~Big_number() = default;

private:
    std::vector<unsigned char> number;
    /// base can be in range [2, 16] on "standard" computer
    unsigned char base;

    void normalise()
    {
        Big_number& ans = *this;
        auto it_end = ans.number.cend() - 1;
        for (auto digit_r = ans.number.begin(); digit_r < it_end; ++digit_r)
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

    [[nodiscard]] Big_number divide(const Big_number& divisor, const bool if_divide) const
    {
        /// dividend / divisor = quotient;
        /// dividend % divisor = remainder;

        Big_number dividend = *this;
        assert (static_cast<bool>(divisor));
        if (dividend < divisor)
        {
            if (if_divide)
                return Big_number();
            else
                return *this;
        }

        Big_number quotient;
        quotient.number.clear();
        Big_number current_part_of_dividend;
        size_t i = dividend.number.size();
        auto p = dividend.number.rbegin();
        do
        {
            current_part_of_dividend = Big_number(*p) + (current_part_of_dividend * 10);
            ++p;
            --i;
        } while (!(current_part_of_dividend >= divisor));
        int quot = current_part_of_dividend.divide_simple(divisor);
        Big_number divisor_on_quot = divisor * quot;
        dividend = dividend - divisor_on_quot * Big_number(10).pow(i);
        quotient.number.push_back(quot);

        bool if_zero = false;
        size_t offset_from_position = divisor_on_quot.number.size();
        for (auto div = divisor_on_quot.number.crbegin(), rem = current_part_of_dividend.number.crbegin();
             *div == *rem && offset_from_position != 1; ++div, ++rem)
        {
            --offset_from_position;
        }
        while (dividend >= divisor || divisor_on_quot == current_part_of_dividend)
        {
            if (divisor_on_quot == current_part_of_dividend)
            {
                dividend.number.push_back(0);
                if_zero = true;
                auto index_zero = this->number.rbegin() + (this->number.size() - i);
                while (index_zero < this->number.rend() && *index_zero == 0)
                {
                    quotient.number.push_back(0);
                    ++index_zero;
                    --i;
                }
                if (index_zero == this->number.rend())
                    goto end;
            }
            p = dividend.number.rbegin();
            current_part_of_dividend = Big_number(*p);

            while (current_part_of_dividend < divisor)
            {
                ++p;
                if (p == dividend.number.rend())// Impossible situation ?
                {
                    assert(false);
                }
                current_part_of_dividend = Big_number(*p) + (current_part_of_dividend * 10);
                if (offset_from_position)
                {
                    if (offset_from_position == 1)
                        --i;
                    --offset_from_position;
                }
                else
                {
                    quotient.number.push_back(0);
                    --i;
                }
            }
            if (if_zero)
            {
                dividend.number.pop_back();
                if_zero = false;
            }
            quot = current_part_of_dividend.divide_simple(divisor);
            divisor_on_quot = divisor * quot;
            assert(!offset_from_position);
            offset_from_position = divisor_on_quot.number.size();/// + 1 -- because it counts jumps - 1 -- because of initialisation
            for (auto div = divisor_on_quot.number.crbegin(), rem = current_part_of_dividend.number.crbegin();
                 *div == *rem && offset_from_position != 1; ++div, ++rem)
            {
                --offset_from_position;
            }
            dividend = dividend - (divisor_on_quot * Big_number(10).pow(i));
            quotient.number.push_back(quot);
        }
        end:
        std::reverse(quotient.number.begin(), quotient.number.end());
        if (if_divide)
            return quotient;
        ///dividend is the remainder
        return dividend;
    }

    [[nodiscard]] int divide_simple(const Big_number& divisor) const
    {
        assert (*this >= divisor);
        int coef = 2;
        Big_number coef_finder = divisor + divisor;
        while (*this >= coef_finder)
        {
            coef_finder += divisor;
            ++coef;
        }
        assert(coef < base);
        return coef - 1;
    }

};

std::ostream& operator<<(std::ostream& ostr, const Big_number& num)
{
    num.view(ostr);
    return ostr;
}

template<>
struct std::is_integral<Big_number>: public true_type{};
#endif //BIG_INT_NUMBERS

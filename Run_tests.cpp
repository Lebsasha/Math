// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_LOG_LEVEL all

//#define NDEBUG 1

#include <boost/test/unit_test.hpp>
#include "Lib/Big_int_numbers.h"
#include "1_Lab/main.cpp"
#include "2_Lab/main.cpp"
#include "3_Lab/main.cpp"
#include "4_Lab/main.cpp"
#include "5_Lab/main.cpp"
#include "6_Lab/main.cpp"
#include "7_Lab/main.cpp"
#include <cfloat>

using namespace std;

BOOST_AUTO_TEST_SUITE(For_Different)

    BOOST_AUTO_TEST_CASE(If_simple_test)
    {
        BOOST_CHECK(if_simple(3) == 1);
        BOOST_CHECK(if_simple(4) == 0);
    }

    BOOST_AUTO_TEST_CASE(Temperature_test)
    {
        BOOST_CHECK_CLOSE(Celsius_to_Farenheigts(-40), -40, DBL_EPSILON);
        BOOST_CHECK_CLOSE(Farenheigts_to_Celsius(-40), -40, DBL_EPSILON);
    }

    BOOST_AUTO_TEST_CASE(Rand_max_and_view_tests)
    {
        short arr[5] = {1, 2, 3, 4, 5};
        fill_random_values(arr, 1, 5, static_cast<short>(0), static_cast<short>(20));
        View(arr, 1, 5);
        BOOST_CHECK(*min_element(arr, arr + 5) == Find_min_by_abs(arr));
        BOOST_CHECK(*max_element(arr, arr + 5) == Find_max_by_abs(arr));
        double arr_2[5] = {1, 2, 3, 4, 5};
        fill_random_values(arr_2, 1, 5, static_cast<double>(0), static_cast<double>(20));
        View(arr_2, 1, 5);
        vector<unsigned long long> vect{1, 2, 3, 4, 5};
        View(vect);
    }

    BOOST_AUTO_TEST_CASE(Pow_test)
    {
        BOOST_CHECK(Pow(4, 5) == 1024);
        BOOST_CHECK(Pow(3, 3) == 27);
    }

BOOST_AUTO_TEST_SUITE_END()

struct Fixture_for_tests
{
    Big_number three;
    Big_number four;

    Fixture_for_tests() : three(3), four(4)
    {}
//    virtual void test_method() final;
};
BOOST_AUTO_TEST_SUITE(For_big_numbers)

    BOOST_AUTO_TEST_CASE(Intercnagable_with_int)
    {
        Big_number num(3);
        int i = num.get_number<int>();
        cout << i << endl;
        BOOST_CHECK(num == Big_number(i));
    }

    BOOST_FIXTURE_TEST_CASE (Operators, Fixture_for_tests)
    {
        BOOST_CHECK((three + four).get_number<int>() == 7);
        BOOST_CHECK((three + four + four + four + three).get_number<int>() == 7 + 7 + 4);
        BOOST_CHECK((four - three) == Big_number(1));
        BOOST_CHECK((Big_number(1006) - Big_number(6)).get_number<int>() == 1000);
        BOOST_CHECK((Big_number(1001) - Big_number(2)).get_number<int>() == 999);
        BOOST_CHECK((four + four + four - three) == Big_number(9));
        BOOST_CHECK((four * three + three * four + four - three * three * three) == Big_number(4 * 3 + 3 * 4 + 4 - 3 * 3 * 3));
        BOOST_CHECK(!static_cast<bool>(four * three + three * four + four - three * three * three - Big_number(1)));
        BOOST_CHECK((four * three).get_number<int>() == 12);
        BOOST_CHECK((Big_number(UINT_MAX) + Big_number(UINT_MAX)) == Big_number(UINT_MAX) * Big_number(2));
        BOOST_CHECK((Big_number(ULLONG_MAX) + Big_number(ULLONG_MAX) + Big_number(ULLONG_MAX)) == Big_number(ULLONG_MAX) * Big_number(3));
        BOOST_CHECK(four.pow(0).get_number<int>() == 1);
        BOOST_CHECK(four.pow(1).get_number<int>() == 4);
        BOOST_CHECK(four.pow(2).get_number<int>() == 16);
        BOOST_CHECK(four.pow(3).get_number<int>() == 64);
        BOOST_CHECK(four == three + Big_number(1));
        BOOST_CHECK(four != three + Big_number(1000));
        BOOST_CHECK(four < three + Big_number(1000));
        BOOST_CHECK(!(four < three));
        BOOST_CHECK(!(four < four));
        BOOST_CHECK(four <= three + Big_number(1));
        BOOST_CHECK(!(four > four + Big_number(10)));
        BOOST_CHECK(four + Big_number(10) > four);
        BOOST_CHECK(!(four > four));
        BOOST_CHECK(four >= four);
        BOOST_CHECK(Big_number(197) / Big_number(15) == Big_number(13));
        BOOST_CHECK(Big_number(1575) / Big_number(15) == Big_number(105));
        BOOST_CHECK(Big_number(150075) / Big_number(15) == Big_number(10005));
        BOOST_CHECK(Big_number(150000) / Big_number(15) == Big_number(10000));
        BOOST_CHECK(four / four == Big_number(1));
        BOOST_CHECK(three / four == Big_number(0));
        BOOST_CHECK(four / three == Big_number(1));
        BOOST_CHECK(Big_number(109007) / Big_number(15) == Big_number(7267));
        Big_number four_m = std::move(four);
        BOOST_CHECK(four_m.get_number<char>() == 4);
        four_m = three;
        BOOST_CHECK(four_m == three);
        BOOST_CHECK(four_m.get_number<unsigned char>() == 3);
        Big_number num(12345);
        num.set_base(16);
        BOOST_CHECK(num.get_base() == 16);
        BOOST_CHECK(num.get_number<int>() == 12345);
        num.set_base(10);
        BOOST_CHECK(num.get_number<int>() == 12345);
        num.view();
        cout << num << endl;
    }

BOOST_AUTO_TEST_SUITE_END()

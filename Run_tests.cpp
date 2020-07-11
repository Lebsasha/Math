//
// Created by alexander on 08/07/2020.
//
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
//#define BOOST_TEST_LOG_LEVEL all
std::string Path = "../Logs/";
//#include "1_Lab/main.cpp"
//#include "2_Lab/main.cpp"
//#include "3_Lab/main.cpp"
//#include "4_Lab/main.cpp"
//#include "5_Lab/main.cpp"
//#include "6_Lab/main.cpp"
//#include "7_Lab/main.cpp"
#include "Different.h"
#include <cfloat>
using namespace std;
//TODO Поддержка существования пути Path
//if(!oFile)
//cout<<error;
BOOST_AUTO_TEST_SUITE(General)
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
    BOOST_AUTO_TEST_CASE(Rand_max_and_view_arr_test)
    {
        short arr[5] = {1, 2, 3, 4, 5};
        fill_random_values(arr, 1, 5, static_cast<short>(0), static_cast<short>(20));
        View(arr, 1, 5);
        BOOST_CHECK(*min_element(arr, arr + 5)==Find_min_by_abs(arr));
        BOOST_CHECK(*max_element(arr, arr + 5)==Find_max_by_abs(arr));
        double arr2[5] = {1, 2, 3, 4, 5};
        fill_random_values(arr2, 1, 5, static_cast<double>(0), static_cast<double>(20));
        View(arr2, 1, 5);
    }
    BOOST_AUTO_TEST_CASE(Vect_view_test)
    {
        vector<unsigned long long> vect {1, 2, 3, 4, 5};
        View(vect);
    }
    BOOST_AUTO_TEST_CASE(Pow_test)
    {
        BOOST_CHECK(Pow(4, 5)==1024);
        BOOST_CHECK(Pow(3, 3)==27);
    }
    BOOST_AUTO_TEST_CASE(Rand_num_test)
    {
        cout<<my_rand_number()<<endl;
    }
    BOOST_AUTO_TEST_CASE(Pause_test)
    {
        Get_Pause();
    }
BOOST_AUTO_TEST_SUITE_END()

//int main (int argc, char** argv)
//{
//}
/**
  *
  * \param
  * \param
  * \return
 * @brief
 * Макросы отличаются уровнем предупреждения (CHECK — ошибка, REQUIRE — критическая ошибка, WARN — предупреждение). Среди макросов есть следующие:

    BOOST_CHECK(условие) — сообщает об ошибке, если условие ложно;
    BOOST_REQUIRE_EQUAL(аргумент_1, аргумент_2) — сообщает о критической ошибке, если аргумент_1 не равен аргумент_2;
    BOOST_WARN_MESSAGE(условие, сообщение) — выводит предупреждение с текстом сообщения, если условие ложно;
    BOOST_CHECK_NO_THROW(выражение) — сообщает об ошибке, если при вычислении выражения вырабатывается исключение;
    BOOST_CHECK_THROW(выражение, исключение) — сообщает об ошибке, если при вычислении выражения не вырабатывается исключение требуемого типа;
    BOOST_CHECK_CLOSE_FRACTION(аргумент_1, аргумент_2, погрешность) — проваливает тест если аргумент_1 не равен аргумент_2 с заданной погрешностью.
    C++ 11 - 201103L
    C++ 11 - 201402L
    C++ 11 - 201703L
    C++ 20 - 201709L
    all 	- report all log messages including the passed test notification
success 	- the same as all
test_suite 	- show test suite messages
message 	- show user messages
warning 	- report warnings issued by user
error 	- report all error conditions
cpp_exception 	- report uncaught c++ exception
system_error 	- report system originated non-fatal errors (for example, timeout or floating point exception)
fatal_error 	- report only user or system originated fatal errors (for example, memory access violation)
nothing 	- does not report any information
 */

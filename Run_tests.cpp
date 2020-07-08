//
// Created by alexander on 08/07/2020.
//
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "1_Lab/main.cpp"
#include "2_Lab/main.cpp"
#include "3_Lab/main.cpp"
#include "4_Lab/main.cpp"
#include "5_Lab/main.cpp"
#include "6_Lab/main.cpp"
#include "7_Lab/main.cpp"
using namespace std;

BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_1)
    BOOST_AUTO_TEST_CASE(Case_for_lab_1)
    {

    }
BOOST_AUTO_TEST_SUITE_END()

//int main (int argc, char** argv)
//{
//
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
 */

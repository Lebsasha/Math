//
// Created by alexander on 08/07/2020.
//
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_LOG_LEVEL all
#include <boost/multiprecision/cpp_int.hpp>

#include "Lib/Big_int_numbers.h"

#define NDEBUG 1

#include <boost/test/unit_test.hpp>
#include "1_Lab/main.cpp"
#include "2_Lab/main.cpp"
#include "3_Lab/main.cpp"
#include "4_Lab/main.cpp"
#include "5_Lab/main.cpp"
#include "6_Lab/main.cpp"
#include "7_Lab/main.cpp"
#include <cfloat>

using namespace std;
//TODO Поддержка существования пути Path
std::string Path = "../Logs/";
//if(!oFile)
//cout<<error;
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
        double arr2[5] = {1, 2, 3, 4, 5};
        fill_random_values(arr2, 1, 5, static_cast<double>(0), static_cast<double>(20));
        View(arr2, 1, 5);
        vector<unsigned long long> vect{1, 2, 3, 4, 5};
        View(vect);
    }

    BOOST_AUTO_TEST_CASE(Pow_test)
    {
        BOOST_CHECK(Pow(4, 5) == 1024);
        BOOST_CHECK(Pow(3, 3) == 27);
    }

    BOOST_AUTO_TEST_CASE(Pause_test)
    {
        Get_Pause();
        int i = 0;
        cin >> i;
        BOOST_CHECK(i != 0);
    }

BOOST_AUTO_TEST_SUITE_END()

struct Fixture_for_tests
{
    Big_number Three;
    Big_number Four;

    Fixture_for_tests() : Three(3), Four(4)
    {}
};
BOOST_AUTO_TEST_SUITE(For_big_numbers)

    BOOST_AUTO_TEST_CASE(Intercnagable_with_int)
    {
        Big_number num(3);
        int i = num.get_Number<int>();
        cout << i << endl;
        BOOST_CHECK(num == Big_number(i));
    }
    //TODO Я пишу это для тренировки
    BOOST_FIXTURE_TEST_CASE (Operators, Fixture_for_tests)
    {
        BOOST_CHECK((Three + Four).get_Number<int>() == 7);
        BOOST_CHECK((Three + Four + Four + Four + Three).get_Number<int>() == 7 + 7 + 4);
        BOOST_CHECK((Four - Three) == Big_number(1));
        BOOST_CHECK((Big_number(1006) - Big_number(6)).get_Number<int>() == 1000);
        BOOST_CHECK((Big_number(1001) - Big_number(2)).get_Number<int>() == 999);
        BOOST_CHECK((Four + Four + Four - Three) == Big_number(9));
        BOOST_CHECK((Four * Three + Three * Four + Four - Three * Three * Three) == Big_number(4 * 3 + 3 * 4 + 4 - 3 * 3 * 3));
        BOOST_CHECK(static_cast<bool>(Four * Three + Three * Four + Four - Three * Three * Three - Big_number(1)) == false);
        BOOST_CHECK((Four * Three).get_Number<int>() == 12);
        BOOST_CHECK((Big_number(UINT_MAX) + Big_number(UINT_MAX)) == Big_number(UINT_MAX) * Big_number(2));
        BOOST_CHECK((Big_number(ULLONG_MAX) + Big_number(ULLONG_MAX) + Big_number(ULLONG_MAX)) == Big_number(ULLONG_MAX) * Big_number(3));
        BOOST_CHECK(Four.pow(0).get_Number<int>() == 1);
        BOOST_CHECK(Four.pow(1).get_Number<int>() == 4);
        BOOST_CHECK(Four.pow(2).get_Number<int>() == 16);
        BOOST_CHECK(Four.pow(3).get_Number<int>() == 64);
        BOOST_CHECK(Four == Three + Big_number(1));
        BOOST_CHECK(Four != Three + Big_number(1000));
        BOOST_CHECK(Four < Three + Big_number(1000));
        BOOST_CHECK(!(Four < Three));
        BOOST_CHECK(!(Four < Four));
        BOOST_CHECK(Four <= Three + Big_number(1));
        BOOST_CHECK(!(Four > Four + Big_number(10)));
        BOOST_CHECK(Four + Big_number(10) > Four);
        BOOST_CHECK(!(Four > Four));
        BOOST_CHECK(Four >= Four);
        BOOST_CHECK(Big_number(197)/Big_number(15) == Big_number(13));
        BOOST_CHECK(Big_number(1575)/Big_number(15) == Big_number(105));
        BOOST_CHECK(Big_number(150075)/Big_number(15) == Big_number(10005));
        BOOST_CHECK(Big_number(150000)/Big_number(15) == Big_number(10000));
        BOOST_CHECK(Four / Four == Big_number(1));
        BOOST_CHECK(Three / Four == Big_number(0));
        BOOST_CHECK(Four / Three == Big_number(1));
        BOOST_CHECK(Big_number(109007)/Big_number(15) == Big_number(7267));
        Big_number Four_m = std::move(Four);
        BOOST_CHECK(Four_m.get_Number<char>() == 4);
        Four_m = Three;
        BOOST_CHECK(Four_m == Three);
        BOOST_CHECK(Four_m.get_Number<unsigned char>() == 3);
        Big_number num(12345);
        num.View();
        num.set_base(16);
        num.View();
        cout << num.get_base() << endl;
        num.set_base(10);
        cout << num << endl;
        BOOST_CHECK(num.get_Number<int>() == 12345);
        num += Big_number(4);
        BOOST_CHECK(num.get_Number<int>() == 12345 + 4);
    }

#include <iterator>

    BOOST_AUTO_TEST_CASE(CCCCASEEEEEEEEE)
    {
        vector<int> a{1, 2, 3};
        auto itr = a.begin();
        std::advance(itr, 2);
        boost::multiprecision::cpp_int num = 5;
        num = num + 6555;

    }

BOOST_AUTO_TEST_SUITE_END()
//int main (int argc, char** argv)
//{
//}
/**
  * std::addressof
  * std::forward
  *  Function 'Divide' should be marked [[nodiscard]]
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

        cout<<__LINE__<<' '<<__FILE__<<' '<<__FUNCTION__<<' '<<__PRETTY_FUNCTION__<<' '<<__cplusplus<<' '<<endl;
        69 /home/alexander/Projects/Num_methods/Run_tests.cpp Fixture_for_tests Fixture_for_tests::Fixture_for_tests() 201703

        std::cout<<std::is_lvalue_reference_v<Big_number&><<std::endl;
        true
        std::cout<<std::is_rvalue_reference_v<Big_number&&><<std::endl;
        true

    template <typename T>
struct SquareMatrix {
 SquareMatrix(std::initializer_list<T> val) v
 : dim{ square_root(val.size()) }, w
 data(dim, std::vector<T>{}) { x
 auto itr = val.begin(); y
 for(size_t row{}; row<dim; row++){
 data[row].assign(itr, itr+dim); z
 itr += dim; {
 }
 }
 T& at(size_t row, size_t col) {
 if (row >= dim || col >= dim)
 throw std::out_of_range{ "Index invalid." }; |
 return data[row][col]; }
 }
 const size_t dim;
private:
 std::vector<std::vector<T>> data;
};









 struct SS
{
    typedef int aaaaa_int;
    typedef int aaaaa_void;
    typedef long aaaaa_long;
};
template<typename SSS>
    typename enable_if<is_void_v<void_t<typename SSS::aaaaa_int, typename SSS::aaaaa_void, typename SSS::aaaaa_long>>, void>::type
    Determine_SS (SSS s)
    {
        cout<<"SS !"<<endl;
    }

template<typename SSS>
typename enable_if<!(is_void_v<void_t<typename SSS::aaaaa_int, typename SSS::aaaaa_void, typename SSS::aaaaa_long>>), void>::type
Determine_SS (SSS s)
{
    cout<<"Not SS"<<endl;
}


1 	. 	Подходит (match) любому символу
2 	^regex 	Регулярное выражение должно совпасть с начального пункта
3 	regex$ 	Регулярное выражение должно совпасть в конце строки.
4 	[abc] 	Создание определения, может подойти к a или b или c.
5 	[abc][vz] 	Создание определения, может подойти к a или b или c, затем следует v или z.
6 	[^abc] 	Когда появляется знак ^ как первый символ в квадратных скобках, он опровергает модель. Это может совпасть с любым символом, кромер a или b или c.
7 	[a-d1-7] 	Рамки: подходит  к строке между a и пунктом d и числом от 1 до 7.
8 	X|Z 	Поиск X или Z.
9 	XZ 	Поиск X а затем Z.
10 	$ 	Проверка завершения строки.
 
11 	\d 	Любое число, краткая форма [0-9]
12 	\D 	Символ неявляющийся числом, краткая форма ​​​​​​​[^0-9]
13 	\s 	Символ пробела, краткая форма​​​​​​​[ \t\n\x0b\r\f]
14 	\S 	Символ неявляющийся пробелом, краткая форма​​​​​​​​​​​​​ [^\s]
15 	\w 	Символ букв, краткая форма​​​​​​​​​​​​​[a-zA-Z_0-9]
16 	\W 	Символ неявляющийся пробелом, краткая форма​​​​​​​​​​​​​[^\w]
17 	\S+ 	Некоторые символы неявляющиеся пробелом (Один или более)
18 	\b 	Символ яляется a-z или A-Z или 0-9 или _, краткая форма​​​​​​​​​​​​​​​​​​​​ [a-zA-Z0-9_].
 
19 	* 	Появляется 0 или много раз, краткая форма​​​​​​​​​​​​​​​​​​​​{0,}
20 	+ 	Появляется 1 или много раз, краткая форма​​​​​​​​​​​​​​​​​​​​​​​​​​​{1,}
21 	? 	Появляется 0 или 1 раз, ? краткая форма​​​​​​​​​​​​​​​​​​​​​​​​​​​​​ {0,1}.
22 	{X} 	Появляется X раз, {}
23 	{X,Y} 	Появляется от X до Y раз.
24 	*? 	* значит появляется 0 или много раз, добавление ? в конце значит поиск самого маленького совпадения.

 */

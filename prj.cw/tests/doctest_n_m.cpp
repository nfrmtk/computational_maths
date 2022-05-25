#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <numerical_methods/n_m.h>
#include "functions.cpp"	
#include <algorithm>

#include <iostream>

void MOut(const MatrixVs& a) {
	std::cout << "  \n";
	for (int i = 0; i < a.rowCount(); i++) {
		for (int j = 0; j < a.colCount(); j++) {
			std::cout << a.at(i, j) << " ";
		}
		std::cout << "\n";
	}
}
void intro() {
    
    std::cout << "1.Matrix algebra presentation.\n";
}

TEST_CASE("matrix algebra") {

    intro();
    std::cout << "1.1 \n Exceptions showcase.\n";
    try {
        MatrixVs a(1, -2);
    }
    catch (const std::exception& e) {
        std::cout << "tried to create matrix with rows amoaunt > 0 and cols amount < 0. caught: ";
        std::cout << e.what() << " \n";
    }
    std::cout << "\n";
    std::cout << "1.2 Matrixes setup\n";
    MatrixVs A(4, 4), B(4, 4);
    for (ptrdiff_t i(0); i < B.rowCount(); i++) {
        for (ptrdiff_t j(0); j < B.colCount(); j++) {
            B(i, j) = i + j + 1 + (i == j) * (i + 1) * (j + 3);
            A(i, j) = (i == j)*2;
        }
    }
    std::cout << "\nmatrix B down below\n";
    MOut(B);
    std::cout << "\nMatrix A down below \n";
    MOut(A);
    std::cout << "\n1.3 Basic arithmetic \n";
    std::cout << "\nA = A + B\nA:";
    A = A + B;
    MOut(A);

    std::cout << "\nA = A - B\nA:";
    A = A - B;
    MOut(A);

    std::cout << "\nA = A * 2\nA:";
    A *= 2;
    MOut(A);

    std::cout << "\nA = A / 3\nA:";
    A /= 3;
    MOut(A);
    std::cout << "\n1.4 Powers of matrixes\n";
    std::cout << "\nB = B^-1\nB:";
    MOut(B.pow(-1));

    std::cout << "\nA = A^-1\nA:";
    MOut(A.pow(-1));

    std::cout << "\nB = B^-1\nB:";
    MOut(B.Inverse());

    std::cout << "\nB = B^3\nB:";
    MOut(B.pow(3));
    std::cout << '\n';
    std::cout << "\n1.5 Cayley_Hamilton theroem showcase\n";
    std::cout << "\nAccording to given theorem f(A)=0, where f is A's characteristic equation\n";
    std::cout << "\nFor given matrix C:\n";
    MatrixVs C(4, 4);
    for (ptrdiff_t i(0); i < C.rowCount(); i++) {
        for (ptrdiff_t j(0); j < C.colCount(); j++) {
            C(i, j) = i - j + 1 + (i == j) * (i + 1) * (j + 3);
        }
    }
    MOut(C);
    std::cout << "\nC's characteristic equation f(x) is: x^4 - 54x^3 + 983x^2 - 7006x + 16020x^0";
    std::cout << "\nlets calculate f(a)\n";
    MOut(MatrixVs(C).pow(4) - 54 * MatrixVs(C).pow(3) + 983 * MatrixVs(C).pow(2) - 7006 * MatrixVs(C) + 16020 * n_m::eye(4));
}


TEST_CASE("limit") {
    std::cout << "2. limit tests\n";
    std::cout << "\n sin(x)/x at 0 limit is ";
	std::cout << n_m::limit(sinlim, 0, n_m::TWO_SIDED) << '\n';
    
    std::cout << "\n x^0.1 at 0 limit is ";
    std::cout << n_m::limit([](double x) {return pow(x, 1.f / 10); }, 0, n_m::TWO_SIDED) << '\n';

    std::cout << "\n e^(1/x) at 0 left limit is ";
    std::cout << n_m::limit([](double x) {return exp(1 / x); }, 0, n_m::LEFT_SIDED) << "\n";

    try {
        double ans(n_m::limit([](double x) {return exp(1 / x); }, 0, n_m::TWO_SIDED));
    }
    catch (const std::logic_error& e) {
        std::cout << e.what() << " <- caught " << typeid(e).name() << "  in limit x->0 of e^(1/x). That is expected behavour, if limit does not exist\n";
    }

}


TEST_CASE("integral") {
    std::cout << "\n 3.Integration test \nn_m::integrate function works best for smooth functions\n and behavours unstably in case of singularities on the [a, b] interval\n";
    std::cout << "sin(x)/x from 0 to 30\n";
    std::cout << n_m::integrate(sinlim, 0, 30);

    std::cout << "\n 1 /(e^x + 1) from 0 to 30\n";
    std::cout << n_m::integrate([](double x) { return 1 / (exp(x) + 1); }, 0, 30);

    std::cout << "\n 1 /(sin(x) + 2) from 0 to 30\n";
    std::cout << n_m::integrate([](double x) { return 1 / (sin(x) + 2); }, 0, 30);

    /*std::cout << "\n exp(x) from 0 to 50\n";
    std::cout << n_m::integrate([](double x) {return exp(x); }, 0, 30);*/

    try {
        std::cout << "\n 1/(sin(x) + 1) from 0 to 30\n";
        std::cout << n_m::integrate([](double x) {return 1 / (1 + sin(x)); }, 0, 30);
    }
    catch (const std::logic_error &e) {
        std::cout << e.what() << "<- exception " << typeid(e).name() << " caught while integrating\n";
    }
    std::cout << "3.1. logarithm test\n";
    try {
        std::cout << "\n ln(x) from 0 to 30\n";
        std::cout << n_m::integrate([](double x) {return log(x); }, 0, 30);
    }
    catch (const std::logic_error &e) {
        std::cout << e.what() << "<- exception " << typeid(e).name() << " caught while integrating\n";
        std::cout << "result is wrong, but if we try integrating from epsilon to 30: \n";
        std::cout << n_m::integrate([](double x) {return log(x); },1000* DBL_EPSILON, 30);
        std::cout << "\nresult is precise\n";   
    }
}
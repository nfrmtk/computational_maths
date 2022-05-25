#pragma once
#ifndef N_M_H_20212912
#define N_M_H_20212912
#include <cstddef>
#include <iostream>
#include <limits>
#include <cmath>
#include <set>
#include <string>


#include <matrixvs/matrixvs.h>
/// side of the limit
typedef short side;
/// small numerical methods namespace. consists of numerical integration and 
/// limit calculating
namespace n_m {
	/// 
	extern double integrate(double (*f)(double), double a, double b);
	extern double limit(double (*f)(double), double x, side SIDE);
	const side RIGHT_SIDED(1);
	const side LEFT_SIDED(-1);
	const side TWO_SIDED(0);
	extern MatrixVs ones(const uint16_t N);
	extern MatrixVs eye(const uint16_t N);

}
#endif //asd
#pragma once
#ifndef MATRIXVS_MATRIXVS_H_20201020
#define MATRIXVS_MATRIXVS_H_20201020

#include <iosfwd>
#include <vector>
#include <cmath>


class MatrixVs {
public:
	MatrixVs() = default;
	MatrixVs(const MatrixVs&) = default;
	MatrixVs(const ptrdiff_t col_count, const ptrdiff_t row_count);
	~MatrixVs() = default;
	MatrixVs& operator=(const MatrixVs&) = default;

	ptrdiff_t rowCount() const noexcept { return n_row_; }
	ptrdiff_t colCount() const noexcept { return n_col_; }


	double& at(const ptrdiff_t i_row, const ptrdiff_t i_col);
	double at(const ptrdiff_t i_row, const ptrdiff_t i_col) const;

	MatrixVs& transpose() noexcept;

	MatrixVs& operator+= (const MatrixVs& lhs);
	MatrixVs& operator-= (const MatrixVs& lhs);

	MatrixVs& operator*= (const double& lhs);
	MatrixVs& operator/= (const double& lhs);

	
	MatrixVs& Inverse();

	double& operator()(ptrdiff_t i, ptrdiff_t j) noexcept;
	double operator()(ptrdiff_t i, ptrdiff_t j) const noexcept;
	double trace();
	double trace(const MatrixVs& lhs);

	double Determinant();

	MatrixVs& multiply(const MatrixVs& rhs);

	MatrixVs& pow(const int16_t Y);

	MatrixVs& gauss_method() noexcept;

	void swapRows(const ptrdiff_t i_first, const ptrdiff_t i_second);
	void swapColumns(const ptrdiff_t i_first, const ptrdiff_t i_second);


private:
	double precision_{std::pow(std::numeric_limits<double>::epsilon(), 1.f/std::min(colCount(),rowCount()))};
	ptrdiff_t n_row_{ 0 };
	ptrdiff_t n_col_{ 0 };
	std::vector<double> data_;
};

inline 	MatrixVs operator+(const MatrixVs& lhs, const MatrixVs& rhs) { return MatrixVs(lhs) += rhs; }
inline 	MatrixVs operator-(const MatrixVs& lhs, const MatrixVs& rhs) { return MatrixVs(lhs) -= rhs; }
inline MatrixVs operator*(const MatrixVs& lhs, double d) { return MatrixVs(lhs) *= d; }
inline MatrixVs operator/(const MatrixVs& lhs, double d) { return MatrixVs(lhs) /= d; }
inline MatrixVs operator*(double d, const MatrixVs& rhs) { return MatrixVs(rhs) *= d; }
inline MatrixVs operator/(double d, const MatrixVs& rhs) { return MatrixVs(rhs) /= d; }

#endif
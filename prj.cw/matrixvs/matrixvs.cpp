#include <matrixvs/matrixvs.h>
#include <stdexcept>
#include <limits>
#include <cmath>
/// constructing empty matrix of given size
/// \exception IE exception if one of argumentss is equal to 0
MatrixVs::MatrixVs(const ptrdiff_t col_count, const ptrdiff_t row_count)
{
	if (col_count <= 0 || row_count <= 0) {
		throw std::invalid_argument("trying to create matrix with non positive rows or cols amount");
	}
	data_ = std::vector<double> (col_count * row_count);
	n_row_ = row_count;
	n_col_ = col_count;

    precision_ = std::pow(std::numeric_limits<double>::epsilon(), 1.f / row_count);
}

/// setter method 
/// \returns pointer to the value
/// \exception IE exception if index is out of matrix
double& MatrixVs::at(const ptrdiff_t i_row, const ptrdiff_t i_col)
{
	if (i_row >= n_row_ || i_col >= n_col_ || i_col < 0 || i_row < 0) {
		throw std::invalid_argument("1 of indexes out of range ");
	}
	double& a = data_[i_row * colCount() + i_col];
	return a;
}

/// getter method
/// \returns const copy of the value
/// \exception throws IE exception if index is out of range
double MatrixVs::at(const ptrdiff_t i_row, const ptrdiff_t i_col) const
{
	if (i_row >= n_row_ || i_col >= n_col_ || i_col < 0 || i_row < 0) {
		throw std::invalid_argument("1 of indexes out of range ");
	}
	const double a(data_[i_row * colCount() + i_col]);
	return a;
}



/// alternative unsafe setter method, UB if i_row or i_col is incorrect.
///
double& MatrixVs::operator()(ptrdiff_t i_row, ptrdiff_t i_col) noexcept
{
	double& a( data_[i_row * colCount() + i_col]);
	return a;
}

/// 
/// alternative unsafe getter method, UB if i_row or i_col is incorrect.
double MatrixVs::operator()(ptrdiff_t i_row, ptrdiff_t i_col) const noexcept
{
	const double a(data_[i_row * colCount() + i_col]);
	return a;
}

/// transposes matrix. 
/// \returns MxN matrix for NxM object
MatrixVs& MatrixVs::transpose() noexcept
{
	MatrixVs newM(n_row_, n_col_);
	for (int i = 0; i < data_.size(); i++) {
		ptrdiff_t row = i / n_col_;
		ptrdiff_t col = i % n_col_;
		newM.at(col, row) = data_[i];
	}
	*this = newM;
	return *this;
}

/// \exception IA if lhs and *this sizes are different
MatrixVs& MatrixVs::operator+=(const MatrixVs& lhs)
{
	if (colCount() != lhs.colCount() || rowCount() != lhs.rowCount()) {
		throw std::invalid_argument("Cant summ matrixes with different sizes.");
	}
	for (ptrdiff_t i{ 0 }; i < colCount(); i++) {
		for (ptrdiff_t j{ 0 }; j < rowCount(); j++) {
			at(i, j) += lhs.at(i, j);
		}
	}
	return *this;
}


/// \exception IA if lhs and *this sizes are different
MatrixVs& MatrixVs::operator-=(const MatrixVs& lhs)
{
	if (colCount() != lhs.colCount() || rowCount() != lhs.rowCount()) {
		throw std::invalid_argument("Cant summ matrixes with different sizes.");
	}
	for (ptrdiff_t i{ 0 }; i < colCount(); i++) {
		for (ptrdiff_t j{ 0 }; j < rowCount(); j++) {
			at(i,j) -= lhs.at(i, j);
		}
	}
	return *this;
}


MatrixVs& MatrixVs::operator*=(const double& lhs)
{
	for (ptrdiff_t i{ 0 }; i < colCount(); i++) {
		for (ptrdiff_t j{ 0 }; j < rowCount(); j++) {
			at(i, j) *= lhs;
		}
	}
	return *this;
}

MatrixVs& MatrixVs::operator/=(const double& lhs)
{
	for (ptrdiff_t i{ 0 }; i < colCount(); i++) {
		for (ptrdiff_t j{ 0 }; j < rowCount(); j++) {
			at(i, j) /= lhs;
		}
	}
	return *this;
}

/// math inverse matrix calculator
/// \exception LE if matrix is not square of if det == 0
MatrixVs& MatrixVs::Inverse()
{
	if (colCount() != rowCount()) {
		throw std::logic_error("Cant compute inverse on non square matrix");
	}
	MatrixVs ans(colCount(), rowCount());
	for (ptrdiff_t i(0); i < rowCount(); i++) {
		for (ptrdiff_t j(0); j < colCount(); j++) {
			ans(i, j) = i == j;
		}
	}
	for (ptrdiff_t subm(0); subm < rowCount(); subm++) {
		ptrdiff_t j(subm);
		for (j = subm; j < colCount(); j++) {
			if (abs(at(subm, j)) > precision_)
			{
				swapColumns(j, subm);
				ans.swapColumns(j, subm);
				j = subm;
				break;
			}
		}
		if (j == colCount()) {
			throw std::logic_error("can't compute inverse of singular matrix");
		}
		for (ptrdiff_t i(subm + 1); i < rowCount(); i++) {
			double coeff(at(i, subm) / at(subm, subm));
			for (ptrdiff_t j(0); j < colCount(); j++) {
				if (j != subm) {
					at(i, j) -= coeff * at(subm, j);
				}
				else {
					at(i, j) = 0;
				}
				ans(i, j) -= coeff * ans(subm, j);
			}
		}
	}
	
	for (ptrdiff_t subm(rowCount() - 1); subm >= 0; subm--) {

		for (ptrdiff_t i(subm - 1); i >= 0; i--) {
			double coeff(at(i, subm) / at(subm, subm));
			at(i, subm) = 0;
			for (ptrdiff_t j(0); j < colCount(); j++) {
				ans(i, j) -= coeff * ans(subm, j);
			}
		}
		for (ptrdiff_t j(0); j < colCount(); j++) {
			ans(subm, j) /= at(subm, subm);
		}
		at(subm, subm) /= at(subm, subm);
	}
	
	*this = ans;
	return *this;
}

/// summ of diagonal elements
/// \exception throws LE if matrix is not square
double MatrixVs::trace()
{
	return trace(*this);
}

/// summ of diagonal elements
/// \exception throws LE if matrix is not square
double MatrixVs::trace(const MatrixVs& lhs)
{
	if (lhs.colCount() != lhs.rowCount()) {
		throw std::logic_error("Cant compute trace on non square matrix");
	}
	double summ = 0;
	for (ptrdiff_t i(0); i < lhs.rowCount(); i++) {
		summ += lhs(i, i);
	}
	return summ;
}

/// determinant of a matrix
/// triangalizes object
/// \exception throws LE if matrix is not square
double  MatrixVs::Determinant()
{
	
	if (colCount() != rowCount()) {
		throw std::logic_error("Cant compute det of non square matrix");
	}
	double det(1);
	MatrixVs preparedMatrix(gauss_method());

	for (ptrdiff_t i(0); i < rowCount(); i++) {
		det *= preparedMatrix(i, i);
	}

	return det;
}


/// right multyplication by matrix rhs
/// \exception throws LE if matrix is not square
MatrixVs& MatrixVs::multiply(const MatrixVs& rhs)
{
	MatrixVs newM(rhs.n_col_, n_row_);
	if (n_col_ != rhs.n_row_) {
		std::string exception = "number of left matrix columns != number of right matrix rows. \n Cant myltiply!";
		throw std::invalid_argument(exception);
	}
	for (int i = 0; i < n_row_; i++) {
		for (int j = 0; j < rhs.n_col_; j++) {
			double sum = 0;
			for (int r = 0; r < n_col_; r++) {
				sum += this->at(i, r) * rhs.at(r, j);
			}
			newM.at(i, j) = sum;
		}
	}


	*this = newM;

	return *this;
}

/// power of a matrix
/// \returns Y'th power of object
/// \exception throws LE exception if Y < 0 and object is singular.
MatrixVs& MatrixVs::pow(const int16_t Y)
{
	if (Y < 0) {
		Inverse();
	}
	if (Y == 0) {
		MatrixVs ans(rowCount(), rowCount());
		for (ptrdiff_t j(0); j < rowCount(); j++) {
			for (ptrdiff_t i(0); i < rowCount(); i++) {
				ans(i, i) = i == j;
			}
		}
		*this = ans;
		return *this;
	}
	MatrixVs cp(*this);
	for (int16_t i(1); i < Y; i++) {
		multiply(cp);
	}
	return *this;
}


/// triangalizes object
/// 
MatrixVs& MatrixVs::gauss_method() noexcept
{
	int mult(1);
	ptrdiff_t zeroRowsAmount(0);
	constexpr double EPS(std::numeric_limits<double>::epsilon());
	MatrixVs &ans(*this);
	for (ptrdiff_t subm(0); subm < rowCount(); subm++) {
		ptrdiff_t j(subm);
		for (j; j < colCount(); j++) {
			if (abs(ans(subm, j)) > EPS) {
				ans.swapColumns(subm, j);
				if ((subm - j) % 2 != 0) mult *= -1;
				break;
			}
		}
		if (j != colCount()) {
			for (ptrdiff_t i(subm + 1); i < ans.rowCount(); i++) {
				double coeff(ans(i, subm) / ans(subm, subm));
				for (j = subm; j < ans.colCount(); j++) {
					if (j != subm) ans(i, j) -= coeff * ans(subm, j);
					else {
						ans(i, j) = 0;
					}
				}
			}
		}
		else
		{
			if (subm < rowCount() - zeroRowsAmount) {		// this line checks if we are not swapping 1 row twice
				subm--;
				zeroRowsAmount++;
				swapRows(subm + 1, rowCount() - zeroRowsAmount - 1);		// swapping last non-zero row with current, if it consist of small enough elements
			}
		}

	}
	for (ptrdiff_t i(0); i < ans.colCount(); i++) {
		ans(0, i) *= mult;
	}
	for (ptrdiff_t i(0); i < ans.rowCount(); i++) {
		for (ptrdiff_t j(0); j < ans.colCount(); j++) {
			if (abs(ans(i, j)) < precision_) {
				ans(i, j) = 0;
			}
		}
	}
	*this = ans;
	return *this;
}


/// swaps i'th and j'th rows of object
///
void MatrixVs::swapRows(const ptrdiff_t i_first, const ptrdiff_t i_second)
{
	if ( i_first < 0 || i_second < 0) {
		throw std::invalid_argument("trying to acess row with negative number");
	}
	if (i_first >= n_row_ || i_second >= n_row_) {
		throw std::invalid_argument("trying to acess row with number > row amount");
	}
	if (i_first != i_second) {
		for (int i = 0; i < n_col_; i++) {
			std::swap(at(i_first, i), at(i_second, i));
		}
	}
}

///  swaps i'th and j'th cols of matrix
///

void MatrixVs::swapColumns(const ptrdiff_t i_first, const ptrdiff_t i_second)
{
	if (i_first >= n_col_ || i_first < 0 || i_second >= n_col_ || i_second < 0) {
		throw std::invalid_argument("trying to acess not existing column");
	}
	if (i_first != i_second) {
		for (int i = 0; i < n_row_; i++) {
			std::swap(at(i, i_first), at(i, i_second));
		}
	}
}




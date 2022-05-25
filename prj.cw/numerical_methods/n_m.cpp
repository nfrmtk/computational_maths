#include <numerical_methods/n_m.h>
namespace n_m {
	namespace {
		double eps = std::numeric_limits<double>::epsilon();
		double simpson_method(double(*f)(double), double a, double b, int n) {
			double reimannSumm = 0;
			double xGap = (b - a) / n;
			for (double i = a + eps; i < b; i+= xGap) {
				reimannSumm += xGap * (f(i) + 4 * f(i + xGap / 2) + f(i + xGap)) / 6;
			}
			return reimannSumm;
		}
		std::set<double> gradient_descent(double(*f) (double), double a, double b, int n) {
			for (a; a < b; a += (b - a) / n) {

			}
		}
		const int32_t MAX_PRECISION(700000);
		const double LIM_DISTANCE (cbrt(std::numeric_limits<double>::epsilon()));
		bool debug(false);
	}

	



	/// numerical integration best suted for integration of smooth function. may behaviour unpredictably dealing with singularities
	/// \exception throws LE if there are singularities on the [a, b] interval
	/// \exception throws LE if ans is too big to be contained in double
	/// \returns signed double
	double integrate(double(*f)(double), double a, double b)
	{
		int mult = 1;
		if (b < a) {
			std::swap(a, b);
			mult = -1;
		}
		double delta = pow(10, -5);
		int n = 11;
		double difference = log(simpson_method(f, a, b, n) / simpson_method(f, a, b, 2 * n));
		
		while (abs(difference) > delta * 15) {
			n *= 2;
			difference = simpson_method(f, a, b, n) - simpson_method(f, a, b, 2 * n);
			if (n > MAX_PRECISION) {
				throw std::logic_error("Integral does not converge ");
			}
		}
				if (!std::isfinite(simpson_method(f, a, b, n))) {
			throw std::logic_error("Answer is too big ");
		}	
		return simpson_method(f, a, b, n) * mult;
	}

	/// function for finding limit. 
	/// \exception throws LE if incorrect side was chosen
	/// \exception throws IA if SIDE is not equal to -1, 0, 1
	/// \exception throws LE if function throws exception at X
	/// \returns signed double
	double limit(double(*f)(double), const double X, const side SIDE)
	{

		bool ifexists(true);
		double ans(std::numeric_limits<double>::infinity());
		try {
			ans = f(X);
		}
		catch ( const std::exception &e) {
			throw std::logic_error("functions is not defined at " + std::to_string(X) + " because of " + e.what());
		}
		if (std::isfinite(ans)) {
			return ans;
		}
		if (SIDE == TWO_SIDED){
			if (abs(f(X+LIM_DISTANCE) - f(X - LIM_DISTANCE) < 10*LIM_DISTANCE)) {
				return (f(X + LIM_DISTANCE) + f(X - LIM_DISTANCE) )/ 2;
			}
			else {
				ifexists = false;
			}
		}
		if (SIDE == LEFT_SIDED) {
			if (abs(f(X - LIM_DISTANCE) - f(X - 2 * LIM_DISTANCE)) < 10*LIM_DISTANCE) {
				return f(X - LIM_DISTANCE);
			}
			else {
				ifexists = false;
			}
		}
		if (SIDE == RIGHT_SIDED) {
			if (abs(f(X + LIM_DISTANCE) - f(X + 2 * LIM_DISTANCE)) < 10*LIM_DISTANCE) {
				return f(X - LIM_DISTANCE);
			}
			else {
				ifexists = false;
			}
		}

		if (!ifexists) {
			throw std::logic_error ("limit doesn't exist");
		}
		throw std::invalid_argument("incorrect side chosen");
	}


	MatrixVs ones(const uint16_t N)
	{
		MatrixVs ones_(N, N);
		for (uint16_t i(0); i < N; i++) {
			for (uint16_t j(0); j < N; j++) {
				ones_(i, j) = 1;
			}
		}
		return ones_;
	}

	MatrixVs eye(const uint16_t N)
	{
		MatrixVs eye_(N, N);
		for (uint16_t i(0); i < N; i++) {
			for (uint16_t j(0); j < N; j++) {
				eye_(i, j) = (i==j);
			}
		}
		return eye_;
	}
}

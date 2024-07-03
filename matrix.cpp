#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/interval.hpp>
#include <random>
#include <iostream>
#include <cmath> // For std::fabs

// Declare namespaces to avoid using namespace directive
namespace ublas = boost::numeric::ublas;
namespace numeric = boost::numeric;

template<typename T>
void gaussianElimination(ublas::matrix<T>& matrix) {
	size_t rows = matrix.size1();
	size_t cols = matrix.size2();

	for (size_t k = 0; k < rows; ++k) { /*
		// Find the maximum element for pivot
		T max = std::fabs(matrix(k,k));
		size_t maxRow = k;
		for (size_t i = k + 1; i < rows; ++i) {
			if (std::fabs(matrix(i, k)) > max) {
				max = std::fabs(matrix(i, k));
				maxRow = i;
			}
		}

		// Swap the maximum row with the current row
		if (k != maxRow)
			for (size_t i = 0; i < cols; ++i)
				std::swap(matrix(k, i), matrix(maxRow, i));*/

		// Perform elimination
		for (size_t i = k + 1; i < rows; ++i) {
			T c = matrix(i, k) / matrix(k, k);
			for (size_t j = k; j < cols; ++j) {
				if (k == j) {
					matrix(i, j) = 0;
				} else {
					matrix(i, j) -= c * matrix(k, j);
				}
			}
		}
	}
}

template<typename T>
ublas::vector<T> backSubstitution(const ublas::matrix<T>& matrix) {
	const size_t n = matrix.size1();
	ublas::vector<T> solution(n);

	for (int i = n - 1; i >= 0; --i) {
		T sum = 0;
		for (size_t j = i + 1; j < n; ++j) {
			sum += matrix(i, j) * solution(j);
		}
		solution(i) = (matrix(i, n) - sum) / matrix(i, i);
		/*
		// Check for division by zero or singular matrix
		if (std::abs(matrix(i, i)) < std::numeric_limits<T>::epsilon()) {
			throw std::runtime_error("Singular matrix or division by zero encountered.");
		}
		*/
	}

	return solution;
}

// Function to generate a square matrix of size 'n' with random float values
ublas::matrix<float> generateRandomMatrix(int n, float range) {
	ublas::matrix<float> matrix(n, n);

	// Random number generation setup
	std::random_device rd;  // Obtain a random number from hardware
	std::mt19937 eng(rd()); // Seed the generator
	std::uniform_real_distribution<> distr(-1 * range, range); // Define the range

	// Fill the matrix with random numbers
	for (size_t i = 0; i < matrix.size1(); ++i) {
		for (size_t j = 0; j < matrix.size2(); ++j) {
			matrix(i, j) = static_cast<float>(distr(eng));
		}
	}

	return matrix;
}

ublas::matrix<numeric::interval<float>> generateRandomMatrixInterval(int n, float range, float error) {
	ublas::matrix<numeric::interval<float>> matrix(n, n);

	// Random number generation setup
	std::random_device rd;  // Obtain a random number from hardware
	std::mt19937 eng(rd()); // Seed the generator
	std::uniform_real_distribution<> distr(-1 * range, range); // Define the range

	// Fill the matrix with random numbers
	for (size_t i = 0; i < matrix.size1(); ++i) {
		for (size_t j = 0; j < matrix.size2(); ++j) {
			float lower = static_cast<float>(distr(eng));
			float upper = lower + error * lower;
			if (lower > upper) std::swap(lower, upper);
			numeric::interval<float> ival(lower, upper);
			matrix(i, j) = ival;
		}
	}

	return matrix;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const boost::numeric::interval<T>& interval) {
	os << '[' << interval.lower() << ", " << interval.upper() << ']';
	return os;
}

template<typename T>
void printMatrix(const ublas::matrix<T>& matrix) {
	for (size_t i = 0; i < matrix.size1(); ++i) {
		for (size_t j = 0; j < matrix.size2(); ++j) {
			std::cout << matrix(i, j) << "\t";
		}
		std::cout << "\n";
	}
}

template<typename T>
void printSolutionVector(const ublas::vector<T>& vector) {
	std::cout << "Solution Vector: [";
	for (size_t i = 0; i < vector.size(); ++i) {
		std::cout << vector(i);
		if (i < vector.size() - 1) std::cout << ", ";
	}
}

template<typename T>
void multiplyMatrixByFactor(ublas::matrix<T>& matrix, float factor) {
	for (size_t i = 0; i < matrix.size1(); ++i) {
		for (size_t j = 0; j < matrix.size2(); ++j) {
			matrix(i, j) *= factor;
		}
	}
}

ublas::vector<float> generateRandomVector(size_t size, float range) {
	ublas::vector<float> vec(size);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1 * range, range);

	for (size_t i = 0; i < size; ++i) {
		vec(i) = static_cast<float>(dis(gen));
	}

	return vec;
}

ublas::vector<numeric::interval<float>> generateRandomVectorInterval(size_t size, float range, float error) {
	ublas::vector<numeric::interval<float>> vec(size);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1 * range, range);

	for (size_t i = 0; i < size; ++i) {
		float lower = static_cast<float>(dis(gen));
		float upper = lower + error * lower;
		if (lower > upper) std::swap(lower, upper);
		numeric::interval<float> ival(lower, upper);
		vec(i) = ival;
	}

	return vec;
}

template<typename T>
void augmentMatrixWithVector(ublas::matrix<T>& matrix, const ublas::vector<T>& vec) {
	size_t rows = matrix.size1();
	size_t cols = matrix.size2();

	// Check if the vector size matches the number of rows in the matrix
	if (vec.size() != rows) {
		throw std::invalid_argument("Vector size must match the number of rows in the matrix.");
	}

	// Resize the matrix to accommodate the vector as an additional column
	matrix.resize(rows, cols + 1, true);

	for (size_t i = 0; i < rows; ++i) {
		matrix(i, cols) = vec(i);
	}
}



int errorFloats() {
	// Generate a random matrix
	ublas::matrix<float> originalMatrix = generateRandomMatrix(100, 0.001);
	ublas::matrix<float> matrix = originalMatrix;

	ublas::vector<float> b = generateRandomVector(100, 0.001);

	// Augment the matrix with the vector b
	augmentMatrixWithVector(matrix, b);

	//std::cout << "\nOriginal Matrix:\n";
	//printMatrix(matrix);

	// Perform Gaussian Elimination
	gaussianElimination(matrix);
	//std::cout << "\nMatrix after Gaussian Elimination:\n";
	//printMatrix(matrix);

	// Perform back substitution
	ublas::vector<float> solution = backSubstitution(matrix);
	std::cout << "\nSolution Vector:\n";
	printSolutionVector(solution);


	matrix = originalMatrix;
	// multiply original matrix by factor
	multiplyMatrixByFactor(matrix, 1.10);

	// Augment the matrix with the vector b
	augmentMatrixWithVector(matrix, b);

	//std::cout << "\n\n\nModified Matrix:\n";
	//printMatrix(matrix);

	// Perform Gaussian Elimination
	gaussianElimination(matrix);
	//std::cout << "\nMatrix after Gaussian Elimination:\n";
	//printMatrix(matrix);

	// Perform back substitution
	ublas::vector<float> modifiedSolution = backSubstitution(matrix);
	std::cout << "\nSolution Vector:\n";
	printSolutionVector(modifiedSolution);

	// Calculate average and maximum error percentages
	float avgError = 0.0;
	float maxError = 0.0;
	for (size_t i = 0; i < solution.size(); ++i) {
		float error = std::fabs((solution(i) - modifiedSolution(i)) / solution(i)) * 100;
		avgError += error;
		if (error > maxError) maxError = error;
	}
	avgError /= solution.size();
	// print results
	std::cout << "\nAverage Error: " << avgError << "%" << std::endl;
	std::cout << "Maximum Error: " << maxError << "%" << std::endl;
	return 0;
}

int errorInterval() {
	// Generate a random matrix
	ublas::matrix<numeric::interval<float>> originalMatrix = generateRandomMatrixInterval(10, 100, 0);
	ublas::matrix<numeric::interval<float>> matrix = originalMatrix;

	ublas::vector<numeric::interval<float>> b = generateRandomVectorInterval(10, 100, 0);

	// Augment the matrix with the vector b
	augmentMatrixWithVector(matrix, b);

	//std::cout << "\nOriginal Matrix:\n";
	//printMatrix(matrix);

	// Perform Gaussian Elimination
	gaussianElimination(matrix);
	//std::cout << "\nMatrix after Gaussian Elimination:\n";
	//printMatrix(matrix);

	// Perform back substitution
	ublas::vector<numeric::interval<float>> solution = backSubstitution(matrix);
	std::cout << "\nSolution Vector:\n";
	printSolutionVector(solution);

	
	// Calculate average and maximum error percentages
	float avgError = 0.0;
	float maxError = 0.0;
	for (size_t i = 0; i < solution.size(); ++i) {
		float error = std::fabs((solution(i).lower() - solution(i).upper()) / solution(i).lower()) * 100;
		avgError += error;
		if (error > maxError) maxError = error;
	}
	avgError /= solution.size();
	// print results
	std::cout << "\nAverage Error: " << avgError << "%" << std::endl;
	std::cout << "Maximum Error: " << maxError << "%" << std::endl;
	return 0;
}

int main() {
	errorInterval();
	return 0;
}
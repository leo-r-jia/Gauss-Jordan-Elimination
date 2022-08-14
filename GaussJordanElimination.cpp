include <iostream> 
#include <fstream>
#include <vector>  
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

class MyMatrix
{
private:

	vector<double> elements;   // array containing data 
	vector<double*> rows;      // array of pointers to start of each row 
	unsigned ncols;            // number of columns 

	void initialise_row_pointers(unsigned nrows);

public:

	MyMatrix(unsigned nrows, unsigned ncols);

	MyMatrix(const MyMatrix& m);

	double get(unsigned irow, unsigned jcol);

	void set(unsigned irow, unsigned jcol, double value);

	void interchange(unsigned irow, unsigned jrow);

	void multiply(unsigned irow, double x);

	void add(unsigned irow, double x, unsigned jrow);

	friend void gaussJordan(MyMatrix& m);

	friend MyMatrix matrixMultiplication(MyMatrix& m1, MyMatrix& m2);

	friend MyMatrix inverse(MyMatrix& A);
};

// MyMatrix constructor that takes matrix dimensions as arguments
MyMatrix::MyMatrix(unsigned nrows, unsigned _ncols)
{
	unsigned nelements = nrows * _ncols;

	ncols = _ncols;

	// resize elements and rows vectors
	elements.resize(nelements);
	rows.resize(nrows);

	initialise_row_pointers(nrows);
}

// MyMatrix copy constructor that takes another matrix as argument
MyMatrix::MyMatrix(const MyMatrix& m)
{
	// check for A=A
	if (this == &m)
		return;

	elements = m.elements;
	ncols = m.ncols;

	rows.resize(m.rows.size());

	initialise_row_pointers(m.rows.size());

}

// function that initialises the row vector
void MyMatrix::initialise_row_pointers(unsigned nrows)
{
	// initialise base pointer
	double* base_pointer = &(elements[0]);

	for (unsigned i = 0; i < nrows; i++)
	{
		unsigned it = ncols * i;
		rows[i] = base_pointer + it;
	}
}

// method returns the element in the ith row and jth column 
double MyMatrix::get(unsigned irow, unsigned jcol)
{
	return *(rows[irow] + jcol);
}

// method sets the element in the ith row and jth column 
void MyMatrix::set(unsigned irow, unsigned jcol, double value)
{
	*(rows[irow] + jcol) = value;
}

// method that performs a row swap
void MyMatrix::interchange(unsigned irow, unsigned jrow)
{
	for (unsigned col = 0; col < ncols; col++)
	{
		// create a copy of value on irow
		double val = get(irow, col);

		// interchange element column by column
		set(irow, col, get(jrow, col));
		set(jrow, col, val);
	}
}

// method that performs a scalar multiplication on a row
void MyMatrix::multiply(unsigned irow, double x)
{
	double val;

	for (unsigned col = 0; col < ncols; col++) {
		// multiply each element column by column
		val = *(rows[irow] + col) * x;
		set(irow, col, val);
	}
}

// method that performs a row sum
void MyMatrix::add(unsigned irow, double x, unsigned jrow)
{
	double val;

	for (unsigned col = 0; col < ncols; col++) {
		// perform add operation column by column
		val = x * get(jrow, col) + get(irow, col);
		set(irow, col, val);
	}
}

// method that performs multiplication of two matrices
MyMatrix matrixMultiplication(MyMatrix& m1, MyMatrix& m2)
{
	// assume square matrices of equal dimensions (for this program)
	unsigned dim = (m1.rows.size());
	assert(dim == m2.rows.size());

	// create product matrix
	MyMatrix productMatrix(dim, dim);

	// calculate matrix multiplication
	for (unsigned i = 0; i < dim; i++)
		for (unsigned j = 0; j < dim; j++)
			for (unsigned k = 0; k < dim; k++)
			{
				double value = productMatrix.get(i, j) + m1.get(i, k) * m2.get(k, j);
				productMatrix.set(i, j, value);
			}

	return productMatrix;
}

// method that performs Gauss-Jordan elimination on a matrix 
void gaussJordan(MyMatrix& m)
{
	unsigned nrows = m.rows.size();
	unsigned ncols = m.ncols;

	// swap rows so largest, leftmost number is on top
	// repeat for every row
	for (unsigned i = 0; i < nrows; i++)
		for (unsigned k = i + 1; k < nrows; k++)
			if (m.get(i, i) < m.get(k, i))
				m.interchange(i, k);

	// perform elementary row operations so that so elements below each row's leading element is zero
	for (unsigned i = 0; i < nrows; i++)
		for (unsigned k = i + 1; k < nrows; k++)
		{
			double x = m.get(k, i) / m.get(i, i);
			m.add(k, -x, i);
		}

	// perform elementary row operations so that elements above each row's leading element is zero
	for (unsigned i = nrows - 1; i > 0; i--)
		for (unsigned k = 0; k < i; k++)
		{
			double x = m.get(k, i) / m.get(i, i);
			m.add(k, -x, i);
		}

	// multiply each row with a scalar to ensure leading element is one
	for (unsigned i = 0; i < nrows; i++)
	{
		double scalar = 1.0 / m.get(i, i);
		m.multiply(i, scalar);
	}
}

// method that finds the inverse of a matrix using an augumented matrix and GaussJordan method 
MyMatrix inverse(MyMatrix& A)
{
	// assuming square matrix
	unsigned nrows = A.rows.size();
	unsigned ncols = A.rows.size() * 2;

	// create an augumented matrix
	MyMatrix augumented(nrows, ncols);

	// assign matrix A to left hand side of augumented matrix
	for (unsigned i = 0; i < A.rows.size(); i++)
		for (unsigned j = 0; j < A.ncols; j++)
			augumented.set(i, j, A.get(i, j));

	// assign identity matrix to right hand side of augumented matrix
	for (unsigned i = 0; i < nrows; i++)
		augumented.set(i, i + nrows, 1.0);

	// perform Gauss-Jordan elimination on the augumented matrix to find inverse(A)
	gaussJordan(augumented);

	// initialize output matrix B
	MyMatrix B(A.rows.size(), A.ncols);

	// set the B matrix to inverse(A)
	for (unsigned i = 0; i < B.rows.size(); i++)
		for (unsigned j = 0; j < B.ncols; j++)
			B.set(i, j, augumented.get(i, j + B.ncols));

	return B;

}

// method that populates matrix with random double precision numbers, computes the inverse, and comp
int runtest(unsigned dim)
{
	// square matrix for this program
	unsigned nrows = dim;
	unsigned ncols = dim;

	MyMatrix matrix(nrows, ncols);

	const double random_max = RAND_MAX;

	// populate matrix with double precision numbers
	for (unsigned i = 0; i < nrows; i++)
		for (unsigned j = 0; j < ncols; j++)
			matrix.set(i, j, rand() / random_max - 0.5);

	// conputes the inverse of matrix
	MyMatrix inv = inverse(matrix);

	// A*A-1 to compute identity matrix
	MyMatrix identity = matrixMultiplication(matrix, inv);

	// set tolerance
	double tolerance = 1e-11;
	unsigned errors = 0;

	// scans identity matrix for precision of values
	for (unsigned i = 0; i < nrows; i++) {
		for (unsigned j = 0; j < ncols; j++) {
			double expect = 0.0;
			if (i == j)
			{
				expect = 1.0;
			}
			// if difference between actual value and expected value is greater than tolerance
			// record as error
			if (fabs(identity.get(i, j) - expect) > tolerance)
			{
				errors++;
			}
		}
	}
	// method returns the number of errors
	return errors;
}

int main()
{
	unsigned ntest = 100000;
	unsigned failures = 0;

	// dimension of square matrix
	unsigned dim = 13;

	// if runtest does not return 0, it has errors
	for (unsigned i = 0; i < ntest; i++)
		if (runtest(dim) != 0)
			failures++;

	cout << failures << " failures out of " << ntest << " tests (" << double(failures) * 100 / ntest << " %)" << endl;

	return 0;
}

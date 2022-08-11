# Gauss-Jordan-Elimination

GaussJordanElimination.cpp implements the Gauss-Jordan elimination algorithm. It has a *MyMatrix* class that holds a vector of elements (which can be visualized as a 2d matrix), a vector of pointers that points to the element that starts each matrix row, and the number of columns.

The method *matrixMultiplication* takes two *MyMatrix* objects as arguments and performs matrix multiplication and returns the product matrix. It assumes two square matrixes of equal dimensions for this algorithm and to keep code as simple as possible. 

The *gaussJordan* method takes one *MyMatrix* object as the argument and performs Gauss-Jordan elimination to reduce it to row-echelon form.  It uses *interchange*, *add*, and m*ultiply* to ensure the matrix satisfies reduced row-echelon form conditions.

The *inverse* method takes one *MyMatrix* object A as the argument, computes its inverse A<sup>-1</sup> and returns it as a new *MyMatrix* object. This is done by setting up an augmented matrix, with the matrix A on the right-hand side and its identity matrix I on the left-hand side.

The *runtest* method takes one unsigned integer value as the dimension that is used to construct square *MyMatrix* “matrices”. It initialises a square matrix A and populates it with random double precision numbers. It calls the inverse method to compute the inverse of the matrix, A<sup>-1</sup>. Then it calls *matrixMultiplication* with A and A<sup>-1</sup> as arguments to compute the identity matrix. Then, the identity matrix is scanned. Each value in the identity matrix is compared against the expected value and a tolerance value. If the difference between the actual value and the expected value is greater than the tolerance, an error counter is incremented by a value of 1. The method returns the number of errors.

Finally, the main method sets the number of *runtest* to call and compute. It counts the number of failures, or no inverse found, when the number returned from *runtest* is greater than zero.

Gauss-Jordan algorithm only utilises row operations, therefore matrices are only partially pivoted. With a tolerance of 1<sup>-11</sup>, approximately 0.9% of all tests are failures.

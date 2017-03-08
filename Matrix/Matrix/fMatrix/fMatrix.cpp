#include "fMatrix.h"
int fMatrix::nMatCount = 0;
// Initinalize constructor.
fMatrix::fMatrix(int n_rows, int n_cols) : rows(n_rows), cols(n_cols)
{
	this->elem = new Float[this->cols*this->rows];
	std::memset(elem, 0, sizeof(Float)*(this->cols*this->rows));
	this->nMatCount++;
}
// Assign constructor.
fMatrix::fMatrix(Float *Array, int n_rows, int n_cols) : rows(n_rows) , cols(n_cols)
{
	//cout << this->rows << n_rows;
	//cout << this->cols << n_cols;
	this->elem = new Float[this->cols*this->rows];
	for (int i = 0; i < this->cols*this->rows; i++) {
		this->elem[i] = Array[i];
	}
	this->nMatCount++;
}
fMatrix::fMatrix(int n_rows, int n_cols, Float *Array) {
	this->cols = n_cols;
	this->rows = n_rows;
	this->elem = new Float[this->cols*this->rows];
	int arrSize = sizeof(Array) / sizeof(Float);
	if (arrSize == (this->cols*this->rows)) {
		for (int i = 0; i < arrSize; i++) {
			this->elem[i] = Array[i];
		}
	}
	else {
		cout << "Syntax Error! \n";
	}
	this->nMatCount++;
}
// Copy constructor.
fMatrix::fMatrix(const fMatrix &A)
{
	fMatrix B(A.rows, A.cols);
	for (int i = 0; i < (B.cols*B.rows); i++) {
		B.elem[i] = A.elem[i];
	}
	this->nMatCount++;
}
// Destructor
fMatrix::~fMatrix(void)
{
	delete[] elem;
	this->nMatCount--;
}
 //1. A+B
fMatrix  operator +  (const fMatrix &A, const fMatrix &B)
{
	fMatrix C(A.rows, A.cols);
	if ((A.rows == B.rows) && (A.cols == B.cols)) {
		
		for (int i = 0; i < (A.rows*A.cols); i++) {
			C.elem [i] = A.elem[i] + B.elem[i];
		}
		
	}
	C.Show();
	return C;
}
// 2. A-B
fMatrix  operator -  (const fMatrix &A)
{	
	fMatrix B(A.rows, A.cols);
	for (int i = 0; i < (A.rows*A.cols); i++) {
		B.elem[i] = -A.elem[i];
	}
	return B;
}
fMatrix  operator -  (const fMatrix &A, const fMatrix &B)
{
	if ((A.rows == B.rows) && (A.cols == B.cols)) {
		fMatrix C(A.rows, A.cols);
		for (int i = 0; i < (A.rows*A.cols); i++) {
			C.elem[i] = A.elem[i] - B.elem[i];
		}
		return C;
	}
	else {
		cout << "The input matrices are invalid\n";
		system("pause");
		exit(1);
	}
}
// 3. c*A or A*c
fMatrix  operator *  (const fMatrix &A, Float B)
{
	fMatrix C(A.rows, A.cols);
	for (int i = 0; i < (A.rows*A.cols); i++) {
		C.elem[i] = A.elem[i] * B;
	}
	return C;
}
fMatrix  operator *  (Float B, const fMatrix &A)
{
	fMatrix C(A.rows, A.cols);
	for (int i = 0; i < (A.rows*A.cols); i++) {
		C.elem[i] = A.elem[i] * B;
	}
	return C;
}
// 4. A/c
fMatrix  operator /  (const fMatrix &A, Float B)
{
	fMatrix C(A.rows, A.cols);
	for (int i = 0; i < (A.rows*A.cols); i++) {
		C.elem[i] = A.elem[i] / B;
	}
	return C;
}
 //5. A*B
fMatrix  operator *  (const fMatrix &A, const fMatrix &B)
{

	return fMatrix();
}
fVector  operator *  (const fMatrix &, const fVector &)
{
	return fVector();
}
fVector  operator *  (const fVector &, const fMatrix &)
{
	return fVector();
}
fMatrix& operator += (fMatrix &, const fMatrix &)
{
	return fMatrix();
}
fMatrix& operator -= (fMatrix &, const fMatrix &)
{
	return fMatrix();
}
fMatrix& operator *= (fMatrix &, Float)
{
	return fMatrix();
}
fMatrix& operator *= (fMatrix &, const fMatrix &)
{
	return fMatrix();
}
fVector& operator *= (fVector &, const fMatrix &)
{
	return fVector();
}
fMatrix& operator /= (fMatrix &, Float)
{
	return fMatrix();
}

/*-------------------------------------------------------------------------*
*                                                                         *
*  FRIEND FUNCTIONS                                                       *
*                                                                         *
*-------------------------------------------------------------------------*/
fMatrix  Transp(const fMatrix &A)	// Transpose of a matrix
{
	fMatrix B(A.rows, A.cols);
	int i = 0;
	for (int col = 0; col < A.cols; col++) {
		for (int row = 0; row < A.rows; row++) {
			int elemA = row * A.cols + col;
			B.elem[i++] = A.elem[elemA];
		}
	}
	return B;
}
//fMatrix  AATransp(const fMatrix &);	// Computes A * Transp(A).
//fMatrix  ATranspA(const fMatrix &);	// Computes Transp(A) * A.
//fMatrix  Outer(const fVector &, const fVector &);  // Computes the outer product of two vectors.
//fMatrix  Identity(int nSize); // Returns an nSizexnSize identity matrix.
//
//fMatrix  Diag(const fVector &);// Returns the square matrix with the elements of the vector d along its diagonal.
//fVector  Diag(const fMatrix &);// Returns the vector consisting of the diagonal elements of the matrix M
//fMatrix  Diag(Float, Float, Float);// Returns the 3 x 3 diagonal matrix with x, y, and z as its diagonal elements.
//
//double   Determinant(const fMatrix &);// Computes the determinant of a square matrix
//double   Trace(const fMatrix &);// Computes the trace of a square matrix
//double   OneNorm(const fMatrix &);// Computes the L1-norm of the matrix A, which is the maximum absolute column sum.
//double   InfNorm(const fMatrix &);// Computes the Inf-norm of the matrix A, which is the maximum absolute row sum.
//
//fMatrix  Inverse(const fMatrix &);// Computes the inverse of a square matrix.
//fMatrix  Cholesky(const fMatrix &);// Computes Cholesky decomposition of a square matrix.	
//fVector  Mean(const fMatrix &);// Computes column mean value of a matrix.	
//fMatrix  Cov(const fMatrix &);// Returns a covariance matrix of a square matrix.
//fMatrix  Cov(const fVector &);// Returns a covariance matrix of a vector, using outer product.
//void     SVDcmp(fMatrix &AU, fVector &W, fMatrix &V); // Computes SVD decomposition of a matrix.
//
//void     ShowMatrix(const fMatrix &);// Print a matrix on screen.
 /*-------------------------------------------------------------------------*
*                                                                         *
*  A S S I G N M E N T    O P E R A T O R S                               *
*                                                                         *
*-------------------------------------------------------------------------*/
// 6. A=B
fMatrix  &fMatrix::operator=(const fMatrix &M)
{
	this->rows = M.rows;
	this->cols = M.cols;
	this->elem = new Float[this->rows*this->cols];

	for (int i = 0; i < (this->cols*this->rows); i++) {
		this->elem[i] = M.elem[i];
	}
	return *this;
}
fMatrix  &fMatrix::operator=(Float s)
{
	this->rows = 1;
	this->cols = 1;
	this->elem = new Float[1];

	this->elem[0] = s;
	return *this;
}

/*-------------------------------------------------------------------------*
*                                                                         *
*  MATRIX OPERATION FUNCTIONS                                             *
*                                                                         *
*-------------------------------------------------------------------------*/
// 7. Swap
//fMatrix &SwapRows(int i1, int i2);
//fMatrix &SwapCols(int j1, int j2);
// 8. Inverse
//fMatrix &Inv(void);
//
//void    SetCol(int col, const fVector &);
//void    SetRow(int row, const fVector &);
//void    SetBlock(int imin, int imax, int jmin, int jmax, const fMatrix &);
//void    SetBlock(int imin, int imax, int jmin, int jmax, const fVector &);
//void    SetSize(int rows, int cols = 0);
//
//fVector  GetCol(int col) const;
//fVector  GetRow(int row) const;
//fMatrix  GetBlock(int imin, int imax, int jmin, int jmax) const;
//
void	fMatrix::Show() const
{
	//cout << "The Matrix is:" << endl;
	for (int row = 0; row < this->rows; row++) {
		for (int col = 0; col < this->cols; col++) {
			int elemA = row*this->cols + col;
			cout << this->elem[elemA] << "\t";
		}
	cout << endl;
	}
}


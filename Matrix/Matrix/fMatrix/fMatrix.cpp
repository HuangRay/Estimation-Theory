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
		//cout << this->elem[i] << "\t";
	}
	cout << endl;
	this->nMatCount++;
}
fMatrix::fMatrix(int n_rows, int n_cols, Float *Array) : rows(n_rows) , cols(n_cols)
{
	//this->cols = n_cols;
	//this->rows = n_rows;
	this->elem = new Float[this->cols*this->rows];
	//int arrSize = sizeof(Array) / sizeof(Float);
	for (int i = 0; i < this->cols*this->rows; i++) {
		this->elem[i] = Array[i];
		//cout << this->elem[i] << "\t";
	}
	this->nMatCount++;
}
// Copy constructor.
fMatrix::fMatrix(const fMatrix &A)
{
	this->rows = A.rows;
	this->cols = A.cols;
	//cout << rows << endl;
	//cout << cols << endl;
	//memccpy(elem, A.elem, rows*cols,);
	this->elem = new Float [this->rows*this->cols];
	for (int i = 0; i < (this->rows * this->cols); i++) {
		this->elem[i] = A.elem[i];
		//cout << this->elem[i] << "\t";
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
	//C.Show();
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
		exit(EXIT_FAILURE);
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
	fMatrix C(A.rows, B.cols);
	if (A.cols == B.rows) {
		for (int i = 0; i < C.rows; i++) {
			for (int j = 0; j < C.cols; j++) {
				C.elem[i*C.cols + j] = A.GetRow(i)*B.GetCol(j);
			}
		}
	}
	else {
		cout << "Input vectors' size are invalid.\n";
	}
	return C;
}
fVector  operator *  (const fMatrix &MatA, const fVector &VecB)
{
	fVector VecC(MatA.rows);
	if (MatA.cols == VecB.Size()) {
		for (int i = 0; i < (VecC.Size()); i++) {
			VecC.Array()[i] = MatA.GetRow(i) * VecB;
		}
	}
	else {
		cout << "Input vectors' size are invalid.\n";
	}
	return VecC;
}
fVector  operator *  (const fVector &VecB, const fMatrix &MatA)
{
	fVector VecC(MatA.cols);
	if (MatA.rows == VecB.Size()) {
		for (int i = 0; i < (VecC.Size()); i++) {
			VecC.Array()[i] = VecB * MatA.GetCol(i);
		}
	}
	else {
		cout << "Input vectors' size are invalid.\n";
	}
	return VecC;
}
fMatrix& operator += (fMatrix &A, const fMatrix &B)
{
	A = A + B;
	return A;
}
fMatrix& operator -= (fMatrix &A, const fMatrix &B)
{
	A = A - B;
	return A;
}
fMatrix& operator *= (fMatrix &A, Float B)
{
	A = A*B;
	return A;
}
fMatrix& operator *= (fMatrix &A, const fMatrix &B)
{
	A = A*B;
	return A;
}
fVector& operator *= (fVector &A, const fMatrix &B)
{
	A = A*B;
	return A;
}
fMatrix& operator /= (fMatrix &A, Float B)
{
	A = A / B;
	return A;
}

/*-------------------------------------------------------------------------*
*                                                                         *
*  FRIEND FUNCTIONS                                                       *
*                                                                         *
*-------------------------------------------------------------------------*/
fMatrix  Transp(const fMatrix &A)	// Transpose of a matrix
{
	fMatrix B(A.cols, A.rows);
	int i = 0;
	for (int col = 0; col < A.cols; col++) {
		for (int row = 0; row < A.rows; row++) {
			int elemA = row * A.cols + col;
			B.elem[i++] = A.elem[elemA];
		}
	}
	return B;
}
fMatrix  AATransp(const fMatrix &A)	// Computes A * Transp(A).
{
	fMatrix B;
	B = A * Transp(A);
	return B;
}
fMatrix  ATranspA(const fMatrix &A)	// Computes Transp(A) * A.
{
	fMatrix B;
	B = Transp(A)* A;
	return B;
}
fMatrix  Outer(const fVector &A, const fVector &B)  // Computes the outer product of two vectors.
{
	fMatrix MatA(A.Array(), A.Size(), 1);
	fMatrix MatB(B.Array(), 1, B.Size());
	MatA *= MatB;
	return MatA;
}
fMatrix  Identity(int nSize) // Returns an nSizexnSize identity matrix.
{
	fMatrix A(nSize, nSize);
	for (int i = 0; i < A.rows*A.cols; i++) {
		int elemA = A.cols*(i / A.rows) + (i / A.rows);
		if (i == elemA) {
			A.elem[i] = 1;
		}
		else {
			A.elem[i] = 0;
		}
	}
	return A;
}
fMatrix  Diag(const fVector &B)// Returns the square matrix with the elements of the vector d along its diagonal.
{
	fMatrix A(B.Size(), B.Size());
	for (int i = 0; i < A.rows*A.cols; i++) {
		int elemA = A.cols*(i / A.cols) + (i / A.cols);
		if (i == elemA) {
			A.elem[i] = B.Array()[(i / A.cols)];
		}
		else {
			A.elem[i] = 0;
		}
	}
	return A;
}
fVector  Diag(const fMatrix &A)// Returns the vector consisting of the diagonal elements of the matrix M
{
	if (A.rows == A.cols) {
		fVector B(A.rows);
		for (int i = 0; i < B.Size(); i++) {
			B.Array()[i] = A.elem[i*A.cols + i];
		}
		return B;
	}
	else {
		cout << "The input matrix is not a square size";
		system("pause");
		exit(EXIT_FAILURE);
	}
		
}
fMatrix  Diag(Float B, Float C, Float D)// Returns the 3 x 3 diagonal matrix with x, y, and z as its diagonal elements.
{
	Float E[3] = { B, C, D };
	fMatrix A(3,3);
	for (int i = 0; i < A.rows*A.cols; i++) {
		int elemA = A.cols*(i / A.cols) + (i / A.cols);
		if (i == elemA) {
			A.elem[i] = E[(i / A.cols)];
		}
		else {
			A.elem[i] = 0;
		}
	}
	return A;
}
double   Determinant(const fMatrix &B)// Computes the determinant of a square matrix
{
	/*
	*			| a1 a2 a3 |		| a5 a6 |		 | a4 a6 |		  | a4 a5 |
	* det(A) =	| a4 a5 a6 | = a1*  | a8 a9 | + a2*  | a7 a9 | + a3*  | a7 a8 |
	*			| a7 a8 a9 |
	*/
	fMatrix A(B);
	double val = 0;
	if (A.cols == A.rows) {
		if (A.cols > 2) {
			fMatrix *arrMat = new fMatrix[A.cols];
			for (int i = 0; i < A.cols; i++) {
				arrMat[i].SetSize(A.rows - 1, A.cols - 1);
			}
			for (int num_matrix = 0; num_matrix < A.cols; num_matrix++) {
				int num_elem = A.cols;
				for (int i = 0; i < arrMat[0].rows*arrMat[0].cols; i++) {
					if (num_matrix == (num_elem%A.cols)) {
						num_elem++;
					}
					arrMat[num_matrix].elem[i] = A.elem[num_elem];
					//arrMat[num_matrix].Show();
					num_elem++;
				}
				//arrMat[num_matrix].Show();
				if (num_matrix % 2 == 1) A.elem[num_matrix] = -A.elem[num_matrix];
				val += A.elem[num_matrix] * Determinant(arrMat[num_matrix]);
			}
			delete[] arrMat;
		}
		else {
			val = A.elem[0] * A.elem[3] - A.elem[1] * A.elem[2];
			//cout << val << endl;
		}
	}
	else {
		cout << "The input matrix is not a square size";
	}
	
	return val;
}
double   Trace(const fMatrix &A)// Computes the trace of a square matrix
{
	if (A.rows == A.cols) {
		fVector B = Diag(A);
		double trace = 0;
		for (int i = 0; i < B.Size(); i++) {
			trace += B.Array()[i];
		}
		return trace;
	}
	else {
		cout << "The input matrix is not a square.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

}
double   OneNorm(const fMatrix &A)// Computes the L1-norm of the matrix A, which is the maximum absolute column sum.
{
	fVector *B = new fVector[A.cols];
	double *onenorm = new double[A.cols];
	for (int i = 0; i < A.cols; i++) {
		B[i] = A.GetCol(i);
		onenorm[i] = OneNorm(B[i]);
	}
	for (int i = 0; i < A.cols; i++) {
		onenorm[0] = (onenorm[i] >= onenorm[0]) ? onenorm[i] : onenorm[0];
	}
	double ans = onenorm[0];
	delete[] B;
	delete[] onenorm;
	return ans;
}
double   InfNorm(const fMatrix &A)// Computes the Inf-norm of the matrix A, which is the maximum absolute row sum.
{
	fVector *B = new fVector[A.rows];
	double *onenorm = new double[A.rows];
	for (int i = 0; i < A.rows; i++) {
		B[i] = A.GetRow(i);
		onenorm[i] = OneNorm(B[i]);
	}
	for (int i = 0; i < A.rows; i++) {
		onenorm[0] = (onenorm[i] >= onenorm[0]) ? onenorm[i] : onenorm[0];
	}
	double ans = onenorm[0];
	delete[] B;
	delete[] onenorm;
	return ans;
}
fMatrix Cofactor(const fMatrix &A)
{
	fMatrix B(A.rows, A.cols);
	fMatrix *cof = new fMatrix[A.rows*A.cols];
	for (int i = 0; i < A.rows*A.cols; i++) {
		cof[i].SetSize(A.rows - 1, A.cols - 1);
		int num_elem = 0;
		for (int j = 0; j < cof[i].rows*cof[i].cols; j++) {
			if ((i/A.cols) == (num_elem /A.cols)) {
				num_elem = num_elem + A.cols;
			}
			if ((i%A.cols) == (num_elem%A.cols)) {
				num_elem = num_elem + 1;
			}
			if ((i / A.cols) == (num_elem / A.cols)) {
				num_elem = num_elem + A.cols;
			}
			cof[i].elem[j] = A.elem[num_elem];
			num_elem++;
		}
		B.elem[i] = pow((-1), i)*Determinant(cof[i]);		
	}
	delete[] cof;
	return B;
}
fMatrix  Inverse(const fMatrix &A)// Computes the inverse of a square matrix.
{
	/*			   
	*			   1
	*	inv(A) = -----  * adj(A), adj(A) = Transpose(Cofactor(A))
	*			 det(A
	*/
	fMatrix B(A);
	B = (1 / Determinant(A)) * Transp(Cofactor(A));
	return B;
}
fMatrix  Cholesky(const fMatrix &A)// Computes Cholesky decomposition of a square matrix.
{
	fMatrix B(A.rows, A.cols);
	double A_elem = 1;
	for (int i = 0; i < A.rows*A.cols; i += A.cols + 1) {
		A_elem = (A.elem[i] < 0) ? A.elem[i] : A_elem;
	}
	if (A_elem > 0) {
		for (int j = 0; j < A.cols; j++) {
			for (int i = 0; i < A.rows; i++) {
				if (i == j) {
					B.elem[i*A.cols + j] = A.elem[i*A.cols + j];
					for (int k = 0; k <= j - 1; k++) {
						B.elem[i*A.cols + j] -= pow(B.elem[j*A.cols + k],2);
					}
					B.elem[i*A.cols + j] = sqrt(B.elem[i*A.cols + j]);
				}
				else if (i > j) {
					B.elem[i*A.cols + j] = A.elem[i*A.cols + j];
					for (int k = 0; k <= j - 1; k++) {
						B.elem[i*A.cols + j] -= B.elem[i*A.cols + k] * B.elem[j*A.cols + k];
					}
					B.elem[i*A.cols + j] = B.elem[i*A.cols + j] / B.elem[j*A.cols + j];
				}			
			}
		}
		return B;
	}
	else {
		cout << "Cholesky ERROR.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
}
fVector  Mean(const fMatrix &A)// Computes column mean value of a matrix.
{
	fVector B(A.cols);
	for (int i = 0; i < B.Size(); i++) {
		B.Array()[i] = Mean(A.GetCol(i));
	}
	return B;
}
fMatrix  Cov(const fMatrix &A)// Returns a covariance matrix of a square matrix.
{
	
	if (A.rows == A.cols) {
		fMatrix B(A.rows, A.cols);
		int n = A.rows;
		//cout << Mean(A.GetCol(0));
		for (int i = 0; i < B.rows; i++) {
			for (int j = 0; j < B.cols; j++) {
				double mean_i = Mean(A.GetCol(i));
				double mean_j = Mean(A.GetCol(j));
				//cout << mean_i << mean_j << endl;
				B.elem[i*B.cols + j] = ((A.GetCol(i) - mean_i)*(A.GetCol(j) - mean_j)) / (n - 1);
			}
		}
		return B;
		
	}
	else {
		cout << "The input matrix's size is invalid\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	
}
fMatrix Cov(const fVector& A)
{
	fVector subMean = A - Mean(A);
	return Outer(subMean, subMean)/(A.Size()-1);
}

//void     SVDcmp(fMatrix &AU, fVector &W, fMatrix &V); // Computes SVD decomposition of a matrix.
//
void     ShowMatrix(const fMatrix &A)// Print a matrix on screen.
{
	A.Show();
}
 /*-------------------------------------------------------------------------*
*                                                                         *
*  A S S I G N M E N T    O P E R A T O R S                               *
*                                                                         *
*-------------------------------------------------------------------------*/
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
fMatrix &fMatrix::SwapRows(int i1, int i2)
{
	fVector A = this->GetRow(i1);
	fVector B = this->GetRow(i2);
	this->SetRow(i1, B);
	this->SetRow(i2, A);
	return *this;
}
fMatrix &fMatrix::SwapCols(int j1, int j2)
{
	fVector A = this->GetCol(j1);
	fVector B = this->GetCol(j2);
	this->SetCol(j1, B);
	this->SetCol(j2, A);
	return *this;
}
fMatrix &fMatrix::Inv(void)
{
	return Inverse(*this);
}
void    fMatrix::SetCol(int col, const fVector &A)
{ 
	if (col < this->cols) {
		if (A.Size() == this->rows) {
			for (int i = 0; i < this->rows; i++) {
				this->elem[i*this->cols + col] = A.Array()[i++];
			}
		}
		else {
			cout << "The input size of the vector is not equal to the  size of the matrix's row\n";
			system("pause");
			exit(EXIT_FAILURE);	
		}
	}
	else {
		cout << "The selection of the column is out of range";
		system("pause");
		exit(EXIT_FAILURE);
	}
	
}
void    fMatrix::SetRow(int row, const fVector &A)
{
	if (row < this->rows) {
		if (A.Size() == this->cols) {
			for (int i = 0; i<this->cols;i++) {
				this->elem[row*this->cols+i] = A.Array()[i];
			}
		}
		else {
			cout << "The input size of the vector is not equal to the  size of the matrix's column\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}
	else {
		cout << "The selection of the row is out of range";
		system("pause");
		exit(EXIT_FAILURE);
	}
}
void    fMatrix::SetBlock(int imin, int imax, int jmin, int jmax, const fMatrix &A)
{
	int row_size = imax - imin + 1;
	int col_size = jmax - jmin + 1;
	if ((row_size > 0) && (col_size > 0) && 
		(row_size <= this->rows) && (col_size <= this->cols) && 
		(row_size == A.rows) && (col_size == A.cols)) {
		int i = 0;
		for (imin; imin <= imax; imin++) {
			for (jmin; jmin <= jmax; jmin++) {
				this->elem[imin*this->cols + jmin] = A.elem[i++];
			}
		}
	}
	else {
		cout << "The input arguments may invalid\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
}

void    fMatrix::SetBlock(int imin, int imax, int jmin, int jmax, const fVector &A)
{
	int row_size = imax - imin + 1;
	int col_size = jmax - jmin + 1;
	if ((row_size > 0) && (col_size > 0) && 
		(row_size <= this->rows) && (col_size <= this->cols) &&
		((row_size == A.Size()) || (col_size == A.Size()))) {
		int i = 0;
		VecType Type = (row_size == 1) ? RowVec : ColVec;
		switch (Type) {
		case RowVec:
			for (jmin; jmin<jmax; jmin++) {
				this->elem[imin*this->cols + jmin] = A.Array()[i++];
			}
			break;
		case ColVec:
			for (imin; imin<imax; imin++) {
				this->elem[imin*this->cols + jmin] = A.Array()[i++];
			}
			break;
		default:
			cout << "It is weird\n";
			break;
		}
	}
	else {
		cout << "The input arguments may invalid\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
}
void    fMatrix::SetSize(int rows, int cols)
{
	this->rows = rows;
	this->cols = cols;
	this->elem = new Float[this->rows*this->cols];
}
fVector  fMatrix::GetCol(int col) const
{
	fVector A(this->rows);
	for (int i = 0; i < this->rows; i++) {
		A.Array()[i] = this->elem[i*this->cols + col];
	}
	return A;
}
fVector  fMatrix::GetRow(int row) const
{
	fVector A(this->cols);
	for (int i = 0; i < this->cols; i++) {
		A.Array()[i] = this->elem[row*this->cols + i];
	}
	return A;
}
fMatrix  fMatrix::GetBlock(int imin, int imax, int jmin, int jmax) const
{
	fMatrix A(imax - imin + 1, jmax - jmin + 1);
	int i = 0;
	for (imin; imin < imax; imin++) {
		for (jmin; jmin < jmax; jmin++) {
			A.elem[i++] = this->elem[imin*this->cols + jmin];
		}
	}
	return A;
}
void	fMatrix::Show() const
{
	//cout << "The Matrix is:" << endl;
	for (int row = 0; row < this->rows; row++) {
		for (int col = 0; col < this->cols; col++) {
			int elemA = row*this->cols + col;
			
			//cout << setw(10) << this->elem[elemA];
			printf("% 10lf", this->elem[elemA]);
		}
	cout << endl;
	}
}

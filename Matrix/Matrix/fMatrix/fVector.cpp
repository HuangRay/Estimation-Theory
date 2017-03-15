#include "fVector.h"

int fVector::nVecCount = 0;
// Initinalize constructor.
fVector::fVector(int size):size(size)
{
	elem = new Float[size];
	//std::cout << size << fVector::size << std::endl;
	std::memset(elem, 0, sizeof(Float)*size);
	//for (int i = 0; i < size; i++) {
	//	std::cout << *(elem + i) << " ";
	//}
	nVecCount++;
	//std::cout << std::endl;
}
// Copy constructor.
fVector::fVector(const fVector &A)
{
	size = A.size;
	elem = new Float[size];
	for (int i = 0; i < size; i++){
		*(elem + i) = *(A.elem+i);
		//std::cout << *(elem + i) << " ";
	}
	nVecCount++;
	//std::cout << std::endl;
}
// Assign constructor.
fVector::fVector(const Float *x, int n)
{
	elem = new Float[n];
	size = n;
	for (int i = 0; i < n; i++) {
		*(elem + i) = *(x + i);
		//std::cout << *(elem + i) << std::endl;
	}
	nVecCount++;
}
fVector::fVector(int n, const Float *x)
{
	elem = new Float[n];
	size = n;
	//std::copy(x, x + n, elem);
	for (int i = 0; i < n; i++) {
		*(elem + i) = *(x + i);
		//std::cout << *(elem + i) << std::endl;
	}
	nVecCount++;
}
fVector::fVector(Float A, Float B)
{
	size = 2;
	elem = new Float(size);
	*(elem + 0) = A;
	*(elem + 1) = B;
	nVecCount++;
}
fVector::fVector(Float A, Float B, Float C)
{
	size = 3;
	elem = new Float(size);
	elem[0] = A;
	elem[1] = B;
	elem[2] = C;
	nVecCount++;
}
fVector::~fVector()
{
	nVecCount--;
	delete [] elem;
}
fVector  operator +  (const fVector &A, const fVector &B)
{
	fVector C(A.size);

	if (A.size == B.size) {
		for (int i = 0; i < A.size; i++) {
			C.elem[i] = A.elem[i] + B.elem[i];
			//std::cout << C.elem[i] << " ";
		}
	}
	else {
		std::cout << "ERROR ! The vectors' size are not equal.";
		exit(1);
	}
	//std::cout << std::endl;
	return C;
}
fVector  operator -  (const fVector &A, const fVector &B) // Binary minus.
{
	fVector C(A.size);
	if (A.size == B.size) {
		for (int i = 0; i < A.size; i++) {
			C.elem[i] = A.elem[i] - B.elem[i];
		}
	}
	else {
		cout << "ERROR ! The vectors' size are not equal.";
		exit(1);
	}
	return C;
}
fVector  operator -  (const fVector &A)  // Unary minus.
{
	fVector B(A.size);
	for (int i = 0; i < B.size; i++) {
		B.elem[i] = -(A.elem[i]);
	}
	return B;
}
fVector  operator -  (const fVector &A, Float B)
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = A.elem[i] - B;
	}
	return C;
}
fVector  operator -  (Float B, const fVector &A) 
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = A.elem[i] - B;
	}
	return C;
}
fVector  operator *  (const fVector &A, Float B)
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = B * A.elem[i];
	}
	return C;
}
fVector  operator *  (Float B, const fVector &A)
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = B * A.elem[i];
	}
	return C;
}
fVector  operator /  (const fVector &A, Float B)
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = A.elem[i]/B;
	}
	return C;
}
fVector  operator /  (const fVector &A, const fVector &B) // Element-wise division
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = A.elem[i] / B.elem[i];
	}
	return C;
}
double   operator *  (const fVector &A, const fVector &B) // Inner-product between two vectors
{
	double C = 0;
	if (A.size == B.size) {
		for (int i = 0; i < A.size; i++) {
			C += A.elem[i] * B.elem[i];
		}
	}
	return C;
}
fVector  operator ^  (const fVector &A, const fVector &B) // Cross-product between two vectors
{
	/*
	* Caculate cross:
	* Cross can only be used in three dimention.
	* A = (a1,a2,a3)	B = (b1,b2,b3)	
	*		| i	 j	 k |	|a2 a3| |a1 a3| |a1 a2|
	* C =	|a1	a2	a3 | = (|b2 b3|,|b1 b3|,|b1 b2| )
	*		|b1	b2	b3 |	
	*/
	fVector C(A.size);
	if ((A.size == 3) && (B.size == 3)) {
		for (int i = 0; i < C.size; i++) {
			C.elem[i] = A.elem[(i + 1) % (A.size)] * B.elem[(i + 2) % (B.size)] -
				A.elem[(i + 2) % (A.size)] * B.elem[(i + 1) % (B.size)];
		}
	}
	else {
		cout << "The input vectors size are invalid\n";
	}
	return C;
}
fVector& operator += (fVector &A, const fVector &B)
{
	A = A + B;
	return A;
}
fVector& operator -= (fVector &A, const fVector &B)
{
	A = A - B;
	return A;
}
fVector& operator *= (fVector &A, Float B)
{
	A = A*B;
	return A;
}
fVector& operator /= (fVector &A, Float B)
{
	A = A / B;
	return A;
}
fVector Min(const fVector &A, const fVector &B)
{
	if(A.size == B.size){
		fVector C(A.size);
		for (int i = 0; i < C.size; i++) {
			C.elem[i] = (A.elem[i] < B.elem[i]) ? A.elem[i] : B.elem[i];
		}
		return C;
	}
	else {
		cout << "The input vectors' size are not equal\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
}
fVector  Max(const fVector &A, const fVector &B) // Element-wise maximum-element extraction between two vectors
{
	if (A.size == B.size) {
		fVector C(A.size);
		for (int i = 0; i < C.size; i++) {
			C.elem[i] = (A.elem[i] > B.elem[i]) ? A.elem[i] : B.elem[i];
		}
		return C;
	}
	else {
		cout << "The input vectors' size are not equal\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
}
double   Dist(const fVector &A, const fVector &B) // Returns two norm distance between two vectors
{
	fVector C;
	if (A.size == B.size) {
		C = A - B;
	}
	else {
		cout << "The vectors' size are not equal";
	}
	return sqrt(C*C);
}
fVector  Normalize(const fVector &A) // Normalizes a vector into an unit vector
{
	return A / sqrt(A*A);
}
double   OneNorm(const fVector &A) // Returns one norm value of a vector
{
	double oneNorm = 0;
	for (int i = 0; i < A.size; i++) {
		oneNorm += fabs(A.elem[i]);
	}
	return oneNorm;
}
double   TwoNorm(const fVector &A) // Returns two norm value of a vector
{
	return sqrt(A*A);
}
double   TwoNormSqr(const fVector &A) // Returns square of the two norm value of a vector
{
	return A*A;
}
fVector  Sqrt(const fVector &A) // Element-wise square root of a vector
{
	fVector B(A);
	for (int i = 0; i < B.size; i++) {
		B.elem[i] = sqrt(B.elem[i]);
	}
	return B;
}
double   Mean(const fVector &A) // Mean value of a vector.
{
	double mean = 0;
	for(int i = 0; i<A.size; i++){
		mean += A.elem[i];
	}
	mean = mean / A.size;
	return mean;
}
double   Var(const fVector &A) // Variance of a vector. 
{
	double mean = Mean(A);
	//cout << mean << endl;
	double var = (A - mean)*(A - mean) / (A.size - 1);
	return var;
}
double   Std(const fVector &A) // Standard derivation of a vector.    	
{
	return sqrt(Var(A));
}
void     ShowVector(const fVector &A, VecType Type)
{
	cout << "The vector is:" << endl;
	if (Type == ColVec) {
		for (int i = 0; i < A.size; i++) {
			cout << A.elem[i] << endl;
		}
	}
	else if (Type == RowVec) {
		for (int i = 0; i < A.size; i++) {
			cout << A.elem[i] << " ";
		}
	}
	else {
		cout << "ERROR!";
	}
}
fVector &fVector::operator=(const fVector &B)
{
	this->size = B.size;
	this->elem = new Float[this->size];
//	SetSize(B.size);
	for (int i = 0; i < size; i++) {
		elem[i] = B.elem[i];
	}
	//cout << elem[0] << " " << elem[1] << " " << elem[2];
	return *this;
}
void    fVector::operator=(Float A)
{
	this->size = 1;
	this->elem = new Float[this->size];
	this->elem[0] = A;
}
void	fVector::SetSize(int size)
{
	fVector::size = size;
	elem = new Float[size];
}

fVector &fVector::Swap(int i,int j)
{
	Float tmp;
	tmp = this->elem[i];
	this->elem[i] = this->elem[j];
	this->elem[j] = tmp;
	return *this;
}
fVector fVector::GetBlock(int i, int j) const
{
	if (j > i) {
		fVector vec(j - i + 1);
		int count = 0;
		for (i; i < j+1; i++) {
			vec.elem[count++] = this->elem[i];
		}
		return vec;
	}
	else {
		cout << "The arguements are ERROR!\n";
		cout << "(const A, const B), B must > A\n";
		return fVector();
	}
}
void    fVector::SetBlock(int i, int j, const fVector &A)
{
	if ((j > i) && (j < this->size) && ((j-i+1) == A.size)) {
		int count = 0;
		for (i; i <= j; i++) {
			this->elem[i] = A.elem[count++];
		}
	}
	else {
		cout << "The arguements are ERROR!\n";
		cout << "(const A, const B), B must be bigger than A\n";
		cout << "B must be smaller than the vector's size \n";
		cout << "The block size must equal to the input vector's size";
	}
}
void	fVector::Show(VecType Type/* = ColVec*/) const
{
	if (Type == ColVec) {
		//cout << "The vector is\n";
		for (int i = 0; i < this->size; i++) {
			cout << elem[i] << endl;
		}
	}
	else if (Type == RowVec) {
		//std::cout << "The vector is\n";
		for (int i = 0; i < this->size; i++) {
			cout << setw(9) << elem[i];
		}
		cout << endl;
	}
	else {
		cout << "ERROR\n";
	}
}



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
		C.elem[i] = B * C.elem[i];
	}
	return C;
}
fVector  operator *  (Float B, const fVector &A)
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = B * C.elem[i];
	}
	return C;
}
fVector  operator /  (const fVector &A, Float B)
{
	fVector C(A.size);
	for (int i = 0; i < C.size; i++) {
		C.elem[i] = C.elem[i]/B;
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
	* Caculate cross
	* A = (a1,a2,a3,...)	B = (b1,b2,b3,...)	
	*		| i	 j	 k |	|a2 a3| |a1 a3| |a1 a2|
	* C =	|a1	a2	a3 | = (|b2 b3|,|b1 b3|,|b1 b2| )
	*		|b1	b2	b3 |	
	*/
	fVector C(A.size);
	if (A.size == B.size) {
		for (int i = 0; i < C.size; i++) {
			//C.elem[i] = ;
		}
	}
	return C;
}
//fVector& operator += (fVector &, const fVector &);
//fVector& operator -= (fVector &, const fVector &);
//fVector& operator *= (fVector &, Float);
//fVector& operator /= (fVector &, Float);
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

//fVector &Swap(int i, int j);
//fVector GetBlock(int i, int j) const;
//void    SetBlock(int i, int j, const fVector &);
void	fVector::Show(VecType Type/* = ColVec*/) const
{
	if (Type == ColVec) {
		cout << "The vector is\n";
		for (int i = 0; i < this->size; i++) {
			cout << elem[i] << endl;
		}
	}
	else if (Type == RowVec) {
		std::cout << "The vector is\n";
		for (int i = 0; i < this->size; i++) {
			cout << elem[i] << " ";
		}
		cout << endl;
	}
	else {
		cout << "ERROR\n";
	}
}
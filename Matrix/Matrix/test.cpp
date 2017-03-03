#include <stdio.h>
#include <stdlib.h>
#include "fMatrix\fVector.h"
#include "Ctest.h"
int main()
{
	Float arrA[3] = { 1,2,3 };
	Float arrB[3] = { 7,8,9 };
	Float constant = 3;
	fVector vA(2);
	fVector vB(3, arrA);
	fVector vC(arrB, 3);
	fVector vD;
	vD = vB + vC;
	vB.Show(ColVec);
	vC.Show(ColVec);
	vD.Show(ColVec);
	(vB - vC).Show(RowVec);
	(-vC).Show(RowVec);
	(vB - constant).Show(RowVec);
	(constant-vB).Show(RowVec);

	Ctest test;
	system("pause");
	return 0;
}
#include "ParamEstimator.h"
CParamEstimator::CParamEstimator()
{

}
CParamEstimator::~CParamEstimator()
{

}
fMatrix* CParamEstimator::SolveOptParam(fVector* pfVecOptParam)
{
	switch (this->EstiMethod){
	case LS:
		return SolveLeastSquares(pfVecOptParam);
		break;
	case WLS:
		return SolveWeightedLeastSquares(pfVecOptParam);
		break;
	case ML:
		return SolveMaximumLikelihood(pfVecOptParam);
		break;
	default:
		std::cout << "Error argument\n";
		system("pause");
		std::exit(EXIT_FAILURE);
		break;
	}
}
fMatrix* CParamEstimator::SolveLeastSquares(fVector* pfVecP_ls)
{
	fVector v;
	fMatrix Vv;
	fMatrix *Var_ls = new fMatrix;
	*pfVecP_ls = Inverse(ATranspA(*this->pLS->pMat_H))*Transp(*this->pLS->pMat_H)*(*this->pLS->pVec_Z);
	v = (*this->pLS->pVec_Z) - (*this->pLS->pMat_H)*(*pfVecP_ls);
	v = v - Mean(v);
	Vv = Outer(v,v) / (v.Size()-1);
	Vv = Diag(Diag(Vv));
	*Var_ls = Inverse(ATranspA(*this->pLS->pMat_H))*Transp(*this->pLS->pMat_H)*Vv*(*this->pLS->pMat_H)*Inverse(ATranspA(*this->pLS->pMat_H));
	this->Mat = *Var_ls;
	delete Var_ls;
	return &this->Mat;
}
fMatrix* CParamEstimator::SolveWeightedLeastSquares(fVector* pfVecP_wls)
{
	fVector v;
	fVector *pfVecP_ls = new fVector;
	fMatrix Vv;
	fMatrix W;
	fMatrix *Var_wls = new fMatrix;
	*pfVecP_ls = Inverse(ATranspA(*this->pWLS->pMat_H))*Transp(*this->pWLS->pMat_H)*(*this->pWLS->pVec_Z);
	v = (*this->pWLS->pVec_Z) - (*this->pWLS->pMat_H)*(*pfVecP_ls);
	W = Diag(Diag(1 / (Identity(Outer(v, v).GetCol(0).Size()) + Outer(v, v))));
	v = v - Mean(v);
	Vv = Outer(v, v) / (v.Size() - 1);
	Vv = Diag(Diag(Vv));
	*pfVecP_wls = Inverse(Transp(*this->pWLS->pMat_H)*W*(*this->pWLS->pMat_H))*(Transp(*this->pWLS->pMat_H))*W*(*this->pWLS->pVec_Z);
	*Var_wls = Inverse(Transp(*this->pWLS->pMat_H)*W*(*this->pWLS->pMat_H))*(Transp(*this->pWLS->pMat_H))*W*Vv*W*(*this->pWLS->pMat_H)*(Inverse(Transp(*this->pWLS->pMat_H)*W*(*this->pWLS->pMat_H))); 
	this->Mat = *Var_wls;
	delete Var_wls;
	delete pfVecP_ls;
	return &this->Mat;
}
fMatrix* CParamEstimator::SolveMaximumLikelihood(fVector* pfVecP_ml)
{
	fVector v;
	fVector *pfVecP_ls = new fVector;
	fMatrix Vv;
	fMatrix W;
	fMatrix *Var_ml = new fMatrix;
	*pfVecP_ls = Inverse(ATranspA(*this->pML->pMat_H))*Transp(*this->pML->pMat_H)*(*this->pML->pVec_Z);
	pfVecP_ls->Show();
	std::cout << endl;
	v = (*this->pML->pVec_Z) - (*this->pML->pMat_H)*(*pfVecP_ls);

	v = v - Mean(v);

	Vv = Outer(v, v) / (v.Size() - 1);

	Vv = Diag(Diag(Vv));

	Vv.Show();
	std::cout << Determinant(Vv) << std::endl;
	//cout << "test" << endl;
	W = Inverse(Vv);
	//W.Show();
	*pfVecP_ml = Inverse(Transp(*this->pML->pMat_H)*W*(*this->pML->pMat_H))*(Transp(*this->pML->pMat_H))*W*(*this->pML->pVec_Z);
	*Var_ml = Inverse(Transp(*this->pML->pMat_H)*Inverse(Vv)*(*this->pML->pMat_H));
	this->Mat = *Var_ml;
	delete Var_ml;
	return &this->Mat;
}
void CParamEstimator::SetParamEstiMethod(ParamEstiMethod Method)
{
	this->EstiMethod = Method;
}
void CParamEstimator::SetMethodParameters(ParamEstiMethod Method, void*	pParam)
{
	switch (Method)
	{
	case LS:
		this->LS_ = *(st_LS_Param*)pParam;
		this->pLS = &this->LS_;
		break;
	case WLS:
		this->WLS_ = *(st_WLS_Param*)pParam;
		this->pWLS = &this->WLS_;
		break;
	case ML:
		this->ML_ = *(st_ML_Param*)pParam;
		this->pML = &this->ML_;
		break;
	default:
		std::cout << "Error argument!\n ";
		break;
	}
}
ParamEstiMethod	CParamEstimator::GetParamEstiMethod(void) const
{
	return this->EstiMethod;
}
void* CParamEstimator::GetMethodParameters(ParamEstiMethod Method) const
{
	switch (Method) {
	case LS:
		return (void*)this->pLS;
		break;
	case WLS:
		return (void*)this->pWLS;
		break;
	case ML:
		return (void*)this->pML;
		break;
	default:
		std::cout << "Error argument!\n ";
		break;
	}
}
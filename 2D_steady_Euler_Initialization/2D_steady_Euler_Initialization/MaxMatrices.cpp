/*
Max4Matrices��������4�������������λ�ô������ֵ���ɵľ���
Max2Matrices��������2�������������λ�ô������ֵ���ɵľ���
*/

#include "Main.h"

using namespace std;
using namespace arma;

arma::mat Max4Matrices(arma::mat A, arma::mat B, arma::mat C, arma::mat D)
{
	umat ZZ11, ZZ12, ZZ21, ZZ22, ZZ31, ZZ32;
	mat MaxMat1, MaxMat2, MaxMat;

	ZZ11 = (A >= B);
	ZZ12 = 1 - ZZ11;
	MaxMat1 = ZZ11 % A + ZZ12 % B;

	ZZ21 = (MaxMat1 >= C);
	ZZ22 = 1 - ZZ21;
	MaxMat2 = ZZ21 % MaxMat1 + ZZ22 % C;

	ZZ31 = (MaxMat2 >= D);
	ZZ32 = 1 - ZZ31;
	MaxMat = ZZ31 % MaxMat2 + ZZ32 % D;

	return MaxMat;
}

arma::mat Max2Matrices(arma::mat A, arma::mat B)
{
	umat ZZ11, ZZ12;
	mat MaxMat;

	ZZ11 = (A >= B);
	ZZ12 = 1 - ZZ11;
	MaxMat = ZZ11 % A + ZZ12 % B;

	return MaxMat;
}
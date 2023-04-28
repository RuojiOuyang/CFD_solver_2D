/*
Convective������������ĵ��ͨ������Q
*/

#include "Main.h"

using namespace std;
using namespace arma;

void Convective()
{
	mat AA1, AA2, BB1, BB2;

	//1.���б߽��ͨ��
	//1.1�߽紦��ͨ��
	//��Boundary_condition�и���

	//1.2�ڲ��߽紦��ͨ��
	//1.2.1�ڲ��߽紦���غ���
	FaceI_U1.cols(1, (long long)Jmax - 2) = 0.5 * (Cell_U1.cols(0, (long long)Jmax - 3) + Cell_U1.cols(1, (long long)Jmax - 2));
	FaceI_U2.cols(1, (long long)Jmax - 2) = 0.5 * (Cell_U2.cols(0, (long long)Jmax - 3) + Cell_U2.cols(1, (long long)Jmax - 2));
	FaceI_U3.cols(1, (long long)Jmax - 2) = 0.5 * (Cell_U3.cols(0, (long long)Jmax - 3) + Cell_U3.cols(1, (long long)Jmax - 2));
	FaceI_U4.cols(1, (long long)Jmax - 2) = 0.5 * (Cell_U4.cols(0, (long long)Jmax - 3) + Cell_U4.cols(1, (long long)Jmax - 2));

	FaceJ_U1.rows(1, (long long)Imax - 2) = 0.5 * (Cell_U1.rows(0, (long long)Imax - 3) + Cell_U1.rows(1, (long long)Imax - 2));
	FaceJ_U2.rows(1, (long long)Imax - 2) = 0.5 * (Cell_U2.rows(0, (long long)Imax - 3) + Cell_U2.rows(1, (long long)Imax - 2));
	FaceJ_U3.rows(1, (long long)Imax - 2) = 0.5 * (Cell_U3.rows(0, (long long)Imax - 3) + Cell_U3.rows(1, (long long)Imax - 2));
	FaceJ_U4.rows(1, (long long)Imax - 2) = 0.5 * (Cell_U4.rows(0, (long long)Imax - 3) + Cell_U4.rows(1, (long long)Imax - 2));

	//1.2.2���ڲ��߽��ϵ��غ��������ڲ��߽��ϵ�ͨ��
	AA1 = FaceI_U1.cols(1, (long long)Jmax - 2) % (arma::pow(FaceI_U2.cols(1, (long long)Jmax - 2) / FaceI_U1.cols(1, (long long)Jmax - 2), 2.0) +
		arma::pow(FaceI_U3.cols(1, (long long)Jmax - 2) / FaceI_U1.cols(1, (long long)Jmax - 2), 2.0)) / 2.0;
	BB1 = (k - 1.0) * (FaceI_U4.cols(1, (long long)Jmax - 2) - AA1);
	FaceI_F1.cols(1, (long long)Jmax - 2) = FaceI_U1.cols(1, (long long)Jmax - 2) % (FaceI_U2.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_nx.cols(1, (long long)Jmax - 2) + FaceI_U3.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_ny.cols(1, (long long)Jmax - 2));
	FaceI_F2.cols(1, (long long)Jmax - 2) = FaceI_U2.cols(1, (long long)Jmax - 2) % (FaceI_U2.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_nx.cols(1, (long long)Jmax - 2) + FaceI_U3.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_ny.cols(1, (long long)Jmax - 2)) + BB1 % FaceI_nx.cols(1, (long long)Jmax - 2);
	FaceI_F3.cols(1, (long long)Jmax - 2) = FaceI_U3.cols(1, (long long)Jmax - 2) % (FaceI_U2.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_nx.cols(1, (long long)Jmax - 2) + FaceI_U3.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_ny.cols(1, (long long)Jmax - 2)) + BB1 % FaceI_ny.cols(1, (long long)Jmax - 2);
	FaceI_F4.cols(1, (long long)Jmax - 2) = (FaceI_U4.cols(1, (long long)Jmax - 2) + BB1) % (FaceI_U2.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_nx.cols(1, (long long)Jmax - 2) + FaceI_U3.cols(1, (long long)Jmax - 2) /
		FaceI_U1.cols(1, (long long)Jmax - 2) % FaceI_ny.cols(1, (long long)Jmax - 2));

	AA2 = FaceJ_U1.rows(1, (long long)Imax - 2) % (arma::pow(FaceJ_U2.rows(1, (long long)Imax - 2) / FaceJ_U1.rows(1, (long long)Imax - 2), 2.0) +
		arma::pow(FaceJ_U3.rows(1, (long long)Imax - 2) / FaceJ_U1.rows(1, (long long)Imax - 2), 2.0)) / 2.0;
	BB2 = (k - 1.0) * (FaceJ_U4.rows(1, (long long)Imax - 2) - AA2);
	FaceJ_F1.rows(1, (long long)Imax - 2) = FaceJ_U1.rows(1, (long long)Imax - 2) % (FaceJ_U2.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_nx.rows(1, (long long)Imax - 2) + FaceJ_U3.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_ny.rows(1, (long long)Imax - 2));
	FaceJ_F2.rows(1, (long long)Imax - 2) = FaceJ_U2.rows(1, (long long)Imax - 2) % (FaceJ_U2.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_nx.rows(1, (long long)Imax - 2) + FaceJ_U3.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_ny.rows(1, (long long)Imax - 2)) + BB2 % FaceJ_nx.rows(1, (long long)Imax - 2);
	FaceJ_F3.rows(1, (long long)Imax - 2) = FaceJ_U3.rows(1, (long long)Imax - 2) % (FaceJ_U2.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_nx.rows(1, (long long)Imax - 2) + FaceJ_U3.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_ny.rows(1, (long long)Imax - 2)) + BB2 % FaceJ_ny.rows(1, (long long)Imax - 2);
	FaceJ_F4.rows(1, (long long)Imax - 2) = (FaceJ_U4.rows(1, (long long)Imax - 2) + BB2) % (FaceJ_U2.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_nx.rows(1, (long long)Imax - 2) + FaceJ_U3.rows(1, (long long)Imax - 2) /
		FaceJ_U1.rows(1, (long long)Imax - 2) % FaceJ_ny.rows(1, (long long)Imax - 2));

	//2.�������б߽��ͨ���󣬿����嵥Ԫ��ͨ����Ҳ����Qij
	Cell_QC1 = FaceI_F1.cols(1, (long long)Jmax - 1) - FaceI_F1.cols(0, (long long)Jmax - 2) + FaceJ_F1.rows(1, (long long)Imax - 1) - FaceJ_F1.rows(0, (long long)Imax - 2);
	Cell_QC2 = FaceI_F2.cols(1, (long long)Jmax - 1) - FaceI_F2.cols(0, (long long)Jmax - 2) + FaceJ_F2.rows(1, (long long)Imax - 1) - FaceJ_F2.rows(0, (long long)Imax - 2);
	Cell_QC3 = FaceI_F3.cols(1, (long long)Jmax - 1) - FaceI_F3.cols(0, (long long)Jmax - 2) + FaceJ_F3.rows(1, (long long)Imax - 1) - FaceJ_F3.rows(0, (long long)Imax - 2);
	Cell_QC4 = FaceI_F4.cols(1, (long long)Jmax - 1) - FaceI_F4.cols(0, (long long)Jmax - 2) + FaceJ_F4.rows(1, (long long)Imax - 1) - FaceJ_F4.rows(0, (long long)Imax - 2);
}
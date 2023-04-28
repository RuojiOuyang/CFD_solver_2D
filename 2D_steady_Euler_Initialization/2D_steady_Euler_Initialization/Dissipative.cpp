/*
Dissipative函数计算人工黏性项，每个中心点处的Dij
Max函数计算四个数中的最大值
*/

#include "Main.h"
#include "MaxMatrices.h"

using namespace std;
using namespace arma;

void Dissipative()
{
	//1.计算lamta
	//1.1边界上的当地通量Jacobi矩阵谱半径的某种近似
	FaceI_lamta.col(0) = Cell_Vol.col(0) / Cell_dt.col(0);
	//cout << Cell_Vol(0, 0) << " " << Cell_dt[numM](0, 0) << endl;
	FaceI_lamta.col((long long)Jmax - 1) = Cell_Vol.col((long long)Jmax - 2) / Cell_dt.col((long long)Jmax - 2);
	FaceJ_lamta.row(0) = 0.5 * (Cell_Vol.row(0) / Cell_dt.row(0) + Cell_Vol.row((long long)Imax - 2) / Cell_dt.row((long long)Imax - 2));
	FaceJ_lamta.row((long long)Imax - 1) = 0.5 * (Cell_Vol.row((long long)Imax - 2) / Cell_dt.row((long long)Imax - 2) + Cell_Vol.row(0) / Cell_dt.row(0));

	//1.2内部边界上的当地通量Jacobi矩阵谱半径的某种近似
	FaceI_lamta.cols(1, (long long)Jmax - 2) = 0.5 * (Cell_Vol.cols(1, (long long)Jmax - 2) / Cell_dt.cols(1, (long long)Jmax - 2) +
		Cell_Vol.cols(0, (long long)Jmax - 3) / Cell_dt.cols(0, (long long)Jmax - 3));
	FaceJ_lamta.rows(1, (long long)Imax - 2) = 0.5 * (Cell_Vol.rows(1, (long long)Imax - 2) / Cell_dt.rows(1, (long long)Imax - 2) +
		Cell_Vol.rows(0, (long long)Imax - 3) / Cell_dt.rows(0, (long long)Imax - 3));

	//2.激波感受因子
	//2.1边界
	//2.1.1下上边界
	muI.col(0) = arma::abs(DummyIS_P.col(0) - 2.0 * DummyIS_P.col(1) + DummyIS_P.col(2)) / arma::abs(DummyIS_P.col(0) + 2.0 * DummyIS_P.col(1) + DummyIS_P.col(2));
	muI.col(1) = arma::abs(Cell_P.col(0) - 2.0 * DummyIS_P.col(0) + DummyIS_P.col(1)) / arma::abs(Cell_P.col(0) + 2.0 * DummyIS_P.col(0) + DummyIS_P.col(1));
	muI.col(2) = arma::abs(Cell_P.col(1) - 2.0 * Cell_P.col(0) + DummyIS_P.col(0)) / arma::abs(Cell_P.col(1) + 2.0 * Cell_P.col(0) + DummyIS_P.col(0));
	muI.col(Jmax) = arma::abs(DummyIN_P.col(0) - 2.0 * Cell_P.col((long long)Jmax - 2) + Cell_P.col((long long)Jmax - 3)) /
		arma::abs(DummyIN_P.col(0) + 2.0 * Cell_P.col((long long)Jmax - 2) + Cell_P.col((long long)Jmax - 3));
	muI.col((long long)Jmax + 1) = arma::abs(DummyIN_P.col(1) - 2.0 * DummyIN_P.col(0) + Cell_P.col((long long)Jmax - 2)) /
		arma::abs(DummyIN_P.col(1) + 2.0 * DummyIN_P.col(0) + Cell_P.col((long long)Jmax - 2));
	muI.col((long long)Jmax + 2) = arma::abs(DummyIN_P.col(2) - 2.0 * DummyIN_P.col(1) + DummyIN_P.col(0)) /
		arma::abs(DummyIN_P.col(2) + 2.0 * DummyIN_P.col(1) + DummyIN_P.col(0));

	//2.1.2左右边界
	muJ.row(0) = arma::abs(DummyJW_P.row(0) - 2.0 * DummyJW_P.row(1) + DummyJW_P.row(2)) / arma::abs(DummyJW_P.row(0) + 2.0 * DummyJW_P.row(1) + DummyJW_P.row(2));
	muJ.row(1) = arma::abs(Cell_P.row(0) - 2.0 * DummyJW_P.row(0) + DummyJW_P.row(1)) / arma::abs(Cell_P.row(0) + 2.0 * DummyJW_P.row(0) + DummyJW_P.row(1));
	muJ.row(2) = arma::abs(Cell_P.row(1) - 2.0 * Cell_P.row(0) + DummyJW_P.row(0)) / arma::abs(Cell_P.row(1) + 2.0 * Cell_P.row(0) + DummyJW_P.row(0));
	muJ.row(Imax) = arma::abs(DummyJE_P.row(0) - 2.0 * Cell_P.row((long long)Imax - 2) + Cell_P.row((long long)Imax - 3)) /
		arma::abs(DummyJE_P.row(0) + 2.0 * Cell_P.row((long long)Imax - 2) + Cell_P.row((long long)Imax - 3));
	muJ.row((long long)Imax + 1) = arma::abs(DummyJE_P.row(1) - 2.0 * DummyJE_P.row(0) + Cell_P.row((long long)Imax - 2)) /
		arma::abs(DummyJE_P.row(1) + 2.0 * DummyJE_P.row(0) + Cell_P.row((long long)Imax - 2));
	muJ.row((long long)Imax + 2) = arma::abs(DummyJE_P.row(2) - 2.0 * DummyJE_P.row(1) + DummyJE_P.row(0)) /
		arma::abs(DummyJE_P.row(2) + 2.0 * DummyJE_P.row(1) + DummyJE_P.row(0));
	//2.2内部
	muI.cols(3, (long long)Jmax - 1) = arma::abs(Cell_P.cols(2, (long long)Jmax - 2) - 2.0 * Cell_P.cols(1, (long long)Jmax - 3) + Cell_P.cols(0, (long long)Jmax - 4)) /
		arma::abs(Cell_P.cols(2, (long long)Jmax - 2) + 2.0 * Cell_P.cols(1, (long long)Jmax - 3) + Cell_P.cols(0, (long long)Jmax - 4));
	muJ.rows(3, (long long)Imax - 1) = arma::abs(Cell_P.rows(2, (long long)Imax - 2) - 2.0 * Cell_P.rows(1, (long long)Imax - 3) + Cell_P.rows(0, (long long)Imax - 4)) /
		arma::abs(Cell_P.rows(2, (long long)Imax - 2) + 2.0 * Cell_P.rows(1, (long long)Imax - 3) + Cell_P.rows(0, (long long)Imax - 4));

	//3.自适应粘性系数ε2和ε4
	double k2, k4;//ε的比例系数
	mat ZEROS1, ZEROS2;//零矩阵备用
	ZEROS1.zeros((long long)Imax - 1, Jmax);
	ZEROS2.zeros(Imax, (long long)Jmax - 1);
	k2 = 0.6;
	k4 = 1.0 / 128.0;
	FaceI_epslion1 = k2 * Max4Matrices(muI.cols(3, (long long)Jmax + 2), muI.cols(2, (long long)Jmax + 1), muI.cols(1, Jmax), muI.cols(0, (long long)Jmax - 1));
	FaceI_epslion2 = Max2Matrices(ZEROS1, k4 - FaceI_epslion1);
	FaceJ_epslion1 = k2 * Max4Matrices(muJ.rows(3, (long long)Imax + 2), muJ.rows(2, (long long)Imax + 1), muJ.rows(1, Imax), muJ.rows(0, (long long)Imax - 1));
	FaceJ_epslion2 = Max2Matrices(ZEROS2, k4 - FaceJ_epslion1);

	//4.计算差分
	//4.1边界
	//4.1.1下上边界
	FaceI_dU1.col(0) = Cell_U1.col(0) - DummyIS_U1.col(0);
	FaceI_dU2.col(0) = Cell_U2.col(0) - DummyIS_U2.col(0);
	FaceI_dU3.col(0) = Cell_U3.col(0) - DummyIS_U3.col(0);
	FaceI_dU4.col(0) = Cell_U4.col(0) - DummyIS_U4.col(0);
	FaceI_dU1.col((long long)Jmax - 1) = DummyIN_U1.col(0) - Cell_U1.col((long long)Jmax - 2);
	FaceI_dU2.col((long long)Jmax - 1) = DummyIN_U2.col(0) - Cell_U2.col((long long)Jmax - 2);
	FaceI_dU3.col((long long)Jmax - 1) = DummyIN_U3.col(0) - Cell_U3.col((long long)Jmax - 2);
	FaceI_dU4.col((long long)Jmax - 1) = DummyIN_U4.col(0) - Cell_U4.col((long long)Jmax - 2);
	FaceI_dddU1.col(0) = Cell_U1.col(1) - 3.0 * Cell_U1.col(0) + 3.0 * DummyIS_U1.col(0) - DummyIS_U1.col(1);
	FaceI_dddU2.col(0) = Cell_U2.col(1) - 3.0 * Cell_U2.col(0) + 3.0 * DummyIS_U2.col(0) - DummyIS_U2.col(1);
	FaceI_dddU3.col(0) = Cell_U3.col(1) - 3.0 * Cell_U3.col(0) + 3.0 * DummyIS_U3.col(0) - DummyIS_U3.col(1);
	FaceI_dddU4.col(0) = Cell_U4.col(1) - 3.0 * Cell_U4.col(0) + 3.0 * DummyIS_U4.col(0) - DummyIS_U4.col(1);
	FaceI_dddU1.col(1) = Cell_U1.col(2) - 3.0 * Cell_U1.col(1) + 3.0 * Cell_U1.col(0) - DummyIS_U1.col(0);
	FaceI_dddU2.col(1) = Cell_U2.col(2) - 3.0 * Cell_U2.col(1) + 3.0 * Cell_U2.col(0) - DummyIS_U2.col(0);
	FaceI_dddU3.col(1) = Cell_U3.col(2) - 3.0 * Cell_U3.col(1) + 3.0 * Cell_U3.col(0) - DummyIS_U3.col(0);
	FaceI_dddU4.col(1) = Cell_U4.col(2) - 3.0 * Cell_U4.col(1) + 3.0 * Cell_U4.col(0) - DummyIS_U4.col(0);
	FaceI_dddU1.col((long long)Jmax - 1) = DummyIN_U1.col(1) - 3.0 * DummyIN_U1.col(0) + 3.0 * Cell_U1.col((long long)Jmax - 2) - Cell_U1.col((long long)Jmax - 3);
	FaceI_dddU2.col((long long)Jmax - 1) = DummyIN_U2.col(1) - 3.0 * DummyIN_U2.col(0) + 3.0 * Cell_U2.col((long long)Jmax - 2) - Cell_U2.col((long long)Jmax - 3);
	FaceI_dddU3.col((long long)Jmax - 1) = DummyIN_U3.col(1) - 3.0 * DummyIN_U3.col(0) + 3.0 * Cell_U3.col((long long)Jmax - 2) - Cell_U3.col((long long)Jmax - 3);
	FaceI_dddU4.col((long long)Jmax - 1) = DummyIN_U4.col(1) - 3.0 * DummyIN_U4.col(0) + 3.0 * Cell_U4.col((long long)Jmax - 2) - Cell_U4.col((long long)Jmax - 3);
	FaceI_dddU1.col((long long)Jmax - 2) = DummyIN_U1.col(0) - 3.0 * Cell_U1.col((long long)Jmax - 2) + 3.0 * Cell_U1.col((long long)Jmax - 3) - Cell_U1.col((long long)Jmax - 4);
	FaceI_dddU2.col((long long)Jmax - 2) = DummyIN_U2.col(0) - 3.0 * Cell_U2.col((long long)Jmax - 2) + 3.0 * Cell_U2.col((long long)Jmax - 3) - Cell_U2.col((long long)Jmax - 4);
	FaceI_dddU3.col((long long)Jmax - 2) = DummyIN_U3.col(0) - 3.0 * Cell_U3.col((long long)Jmax - 2) + 3.0 * Cell_U3.col((long long)Jmax - 3) - Cell_U3.col((long long)Jmax - 4);
	FaceI_dddU4.col((long long)Jmax - 2) = DummyIN_U4.col(0) - 3.0 * Cell_U4.col((long long)Jmax - 2) + 3.0 * Cell_U4.col((long long)Jmax - 3) - Cell_U4.col((long long)Jmax - 4);
	//4.1.2左右边界
	FaceJ_dU1.row(0) = Cell_U1.row(0) - DummyJW_U1.row(0);
	FaceJ_dU2.row(0) = Cell_U2.row(0) - DummyJW_U2.row(0);
	FaceJ_dU3.row(0) = Cell_U3.row(0) - DummyJW_U3.row(0);
	FaceJ_dU4.row(0) = Cell_U4.row(0) - DummyJW_U4.row(0);
	FaceJ_dU1.row((long long)Imax - 1) = DummyJE_U1.row(0) - Cell_U1.row((long long)Imax - 2);
	FaceJ_dU2.row((long long)Imax - 1) = DummyJE_U2.row(0) - Cell_U2.row((long long)Imax - 2);
	FaceJ_dU3.row((long long)Imax - 1) = DummyJE_U3.row(0) - Cell_U3.row((long long)Imax - 2);
	FaceJ_dU4.row((long long)Imax - 1) = DummyJE_U4.row(0) - Cell_U4.row((long long)Imax - 2);
	FaceJ_dddU1.row(0) = Cell_U1.row(1) - 3.0 * Cell_U1.row(0) + 3.0 * DummyJW_U1.row(0) - DummyJW_U1.row(1);
	FaceJ_dddU2.row(0) = Cell_U2.row(1) - 3.0 * Cell_U2.row(0) + 3.0 * DummyJW_U2.row(0) - DummyJW_U2.row(1);
	FaceJ_dddU3.row(0) = Cell_U3.row(1) - 3.0 * Cell_U3.row(0) + 3.0 * DummyJW_U3.row(0) - DummyJW_U3.row(1);
	FaceJ_dddU4.row(0) = Cell_U4.row(1) - 3.0 * Cell_U4.row(0) + 3.0 * DummyJW_U4.row(0) - DummyJW_U4.row(1);
	FaceJ_dddU1.row(1) = Cell_U1.row(2) - 3.0 * Cell_U1.row(1) + 3.0 * Cell_U1.row(0) - DummyJW_U1.row(0);
	FaceJ_dddU2.row(1) = Cell_U2.row(2) - 3.0 * Cell_U2.row(1) + 3.0 * Cell_U2.row(0) - DummyJW_U2.row(0);
	FaceJ_dddU3.row(1) = Cell_U3.row(2) - 3.0 * Cell_U3.row(1) + 3.0 * Cell_U3.row(0) - DummyJW_U3.row(0);
	FaceJ_dddU4.row(1) = Cell_U4.row(2) - 3.0 * Cell_U4.row(1) + 3.0 * Cell_U4.row(0) - DummyJW_U4.row(0);
	FaceJ_dddU1.row((long long)Imax - 1) = DummyJE_U1.row(1) - 3.0 * DummyJE_U1.row(0) + 3 * Cell_U1.row((long long)Imax - 2) - Cell_U1.row((long long)Imax - 3);
	FaceJ_dddU2.row((long long)Imax - 1) = DummyJE_U2.row(1) - 3.0 * DummyJE_U2.row(0) + 3 * Cell_U2.row((long long)Imax - 2) - Cell_U2.row((long long)Imax - 3);
	FaceJ_dddU3.row((long long)Imax - 1) = DummyJE_U3.row(1) - 3.0 * DummyJE_U3.row(0) + 3 * Cell_U3.row((long long)Imax - 2) - Cell_U3.row((long long)Imax - 3);
	FaceJ_dddU4.row((long long)Imax - 1) = DummyJE_U4.row(1) - 3.0 * DummyJE_U4.row(0) + 3 * Cell_U4.row((long long)Imax - 2) - Cell_U4.row((long long)Imax - 3);
	FaceJ_dddU1.row((long long)Imax - 2) = DummyJE_U1.row(0) - 3.0 * Cell_U1.row((long long)Imax - 2) + 3.0 * Cell_U1.row((long long)Imax - 3) - Cell_U1.row((long long)Imax - 4);
	FaceJ_dddU2.row((long long)Imax - 2) = DummyJE_U2.row(0) - 3.0 * Cell_U2.row((long long)Imax - 2) + 3.0 * Cell_U2.row((long long)Imax - 3) - Cell_U2.row((long long)Imax - 4);
	FaceJ_dddU3.row((long long)Imax - 2) = DummyJE_U3.row(0) - 3.0 * Cell_U3.row((long long)Imax - 2) + 3.0 * Cell_U3.row((long long)Imax - 3) - Cell_U3.row((long long)Imax - 4);
	FaceJ_dddU4.row((long long)Imax - 2) = DummyJE_U4.row(0) - 3.0 * Cell_U4.row((long long)Imax - 2) + 3.0 * Cell_U4.row((long long)Imax - 3) - Cell_U4.row((long long)Imax - 4);

	//4.2内部
	FaceI_dU1.cols(1, (long long)Jmax - 2) = Cell_U1.cols(1, (long long)Jmax - 2) - Cell_U1.cols(0, (long long)Jmax - 3);
	FaceI_dU2.cols(1, (long long)Jmax - 2) = Cell_U2.cols(1, (long long)Jmax - 2) - Cell_U2.cols(0, (long long)Jmax - 3);
	FaceI_dU3.cols(1, (long long)Jmax - 2) = Cell_U3.cols(1, (long long)Jmax - 2) - Cell_U3.cols(0, (long long)Jmax - 3);
	FaceI_dU4.cols(1, (long long)Jmax - 2) = Cell_U4.cols(1, (long long)Jmax - 2) - Cell_U4.cols(0, (long long)Jmax - 3);
	FaceJ_dU1.rows(1, (long long)Imax - 2) = Cell_U1.rows(1, (long long)Imax - 2) - Cell_U1.rows(0, (long long)Imax - 3);
	FaceJ_dU2.rows(1, (long long)Imax - 2) = Cell_U2.rows(1, (long long)Imax - 2) - Cell_U2.rows(0, (long long)Imax - 3);
	FaceJ_dU3.rows(1, (long long)Imax - 2) = Cell_U3.rows(1, (long long)Imax - 2) - Cell_U3.rows(0, (long long)Imax - 3);
	FaceJ_dU4.rows(1, (long long)Imax - 2) = Cell_U4.rows(1, (long long)Imax - 2) - Cell_U4.rows(0, (long long)Imax - 3);
	FaceI_dddU1.cols(2, (long long)Jmax - 3) = Cell_U1.cols(3, (long long)Jmax - 2) - 3.0 * Cell_U1.cols(2, (long long)Jmax - 3) +
		3.0 * Cell_U1.cols(1, (long long)Jmax - 4) - Cell_U1.cols(0, (long long)Jmax - 5);
	FaceI_dddU2.cols(2, (long long)Jmax - 3) = Cell_U2.cols(3, (long long)Jmax - 2) - 3.0 * Cell_U2.cols(2, (long long)Jmax - 3) +
		3.0 * Cell_U2.cols(1, (long long)Jmax - 4) - Cell_U2.cols(0, (long long)Jmax - 5);
	FaceI_dddU3.cols(2, (long long)Jmax - 3) = Cell_U3.cols(3, (long long)Jmax - 2) - 3.0 * Cell_U3.cols(2, (long long)Jmax - 3) +
		3.0 * Cell_U3.cols(1, (long long)Jmax - 4) - Cell_U3.cols(0, (long long)Jmax - 5);
	FaceI_dddU4.cols(2, (long long)Jmax - 3) = Cell_U4.cols(3, (long long)Jmax - 2) - 3.0 * Cell_U4.cols(2, (long long)Jmax - 3) +
		3.0 * Cell_U4.cols(1, (long long)Jmax - 4) - Cell_U4.cols(0, (long long)Jmax - 5);
	FaceJ_dddU1.rows(2, (long long)Imax - 3) = Cell_U1.rows(3, (long long)Imax - 2) - 3.0 * Cell_U1.rows(2, (long long)Imax - 3) +
		3.0 * Cell_U1.rows(1, (long long)Imax - 4) - Cell_U1.rows(0, (long long)Imax - 5);
	FaceJ_dddU2.rows(2, (long long)Imax - 3) = Cell_U2.rows(3, (long long)Imax - 2) - 3.0 * Cell_U2.rows(2, (long long)Imax - 3) +
		3.0 * Cell_U2.rows(1, (long long)Imax - 4) - Cell_U2.rows(0, (long long)Imax - 5);
	FaceJ_dddU3.rows(2, (long long)Imax - 3) = Cell_U3.rows(3, (long long)Imax - 2) - 3.0 * Cell_U3.rows(2, (long long)Imax - 3) +
		3.0 * Cell_U3.rows(1, (long long)Imax - 4) - Cell_U3.rows(0, (long long)Imax - 5);
	FaceJ_dddU4.rows(2, (long long)Imax - 3) = Cell_U4.rows(3, (long long)Imax - 2) - 3.0 * Cell_U4.rows(2, (long long)Imax - 3) +
		3.0 * Cell_U4.rows(1, (long long)Imax - 4) - Cell_U4.rows(0, (long long)Imax - 5);

	//5.所有边界上的人工黏性项
	FaceI_D1 = FaceI_lamta % (FaceI_epslion1 % FaceI_dU1 - FaceI_epslion2 % FaceI_dddU1);
	FaceI_D2 = FaceI_lamta % (FaceI_epslion1 % FaceI_dU2 - FaceI_epslion2 % FaceI_dddU2);
	FaceI_D3 = FaceI_lamta % (FaceI_epslion1 % FaceI_dU3 - FaceI_epslion2 % FaceI_dddU3);
	FaceI_D4 = FaceI_lamta % (FaceI_epslion1 % FaceI_dU4 - FaceI_epslion2 % FaceI_dddU4);
	FaceJ_D1 = FaceJ_lamta % (FaceJ_epslion1 % FaceJ_dU1 - FaceJ_epslion2 % FaceJ_dddU1);
	FaceJ_D2 = FaceJ_lamta % (FaceJ_epslion1 % FaceJ_dU2 - FaceJ_epslion2 % FaceJ_dddU2);
	FaceJ_D3 = FaceJ_lamta % (FaceJ_epslion1 % FaceJ_dU3 - FaceJ_epslion2 % FaceJ_dddU3);
	FaceJ_D4 = FaceJ_lamta % (FaceJ_epslion1 % FaceJ_dU4 - FaceJ_epslion2 % FaceJ_dddU4);

	//6.控制体人工黏性
	Cell_D1 = FaceI_D1.cols(1, (long long)Jmax - 1) - FaceI_D1.cols(0, (long long)Jmax - 2) + FaceJ_D1.rows(1, (long long)Imax - 1) - FaceJ_D1.rows(0, (long long)Imax - 2);
	Cell_D2 = FaceI_D2.cols(1, (long long)Jmax - 1) - FaceI_D2.cols(0, (long long)Jmax - 2) + FaceJ_D2.rows(1, (long long)Imax - 1) - FaceJ_D2.rows(0, (long long)Imax - 2);
	Cell_D3 = FaceI_D3.cols(1, (long long)Jmax - 1) - FaceI_D3.cols(0, (long long)Jmax - 2) + FaceJ_D3.rows(1, (long long)Imax - 1) - FaceJ_D3.rows(0, (long long)Imax - 2);
	Cell_D4 = FaceI_D4.cols(1, (long long)Jmax - 1) - FaceI_D4.cols(0, (long long)Jmax - 2) + FaceJ_D4.rows(1, (long long)Imax - 1) - FaceJ_D4.rows(0, (long long)Imax - 2);
}
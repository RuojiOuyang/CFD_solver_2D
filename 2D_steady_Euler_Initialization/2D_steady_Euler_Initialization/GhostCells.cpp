/*
GhostCells函数，处理虚元，为人工黏性项的计算做铺垫
*/

#include "Main.h"

using namespace std;
using namespace arma;

void GhostCells()
{
	int nd;//迭代量，迭代虚元数量时使用
	mat V_f, V_a, V_b;

	//1.算出各个边界处的虚元
	//1.1物面边界的虚元
	DummyIS_Den.col(0) = Cell_Den.col(0);
	DummyIS_E.col(0) = Cell_E.col(0);
	DummyIS_P.col(0) = Cell_P.col(0);
	DummyIS_T.col(0) = Cell_T.col(0);
	V_f = Cell_u.col(0) % FaceI_nx.col(0) + Cell_v.col(0) % FaceI_ny.col(0);
	V_a = V_f % FaceI_nx.col(0) / (FaceI_nx.col(0) % FaceI_nx.col(0) + FaceI_ny.col(0) % FaceI_ny.col(0));
	V_b = V_f % FaceI_ny.col(0) / (FaceI_nx.col(0) % FaceI_nx.col(0) + FaceI_ny.col(0) % FaceI_ny.col(0));
	DummyIS_u.col(0) = Cell_u.col(0) - 2 * V_a;
	DummyIS_v.col(0) = Cell_v.col(0) - 2 * V_b;
	DummyIS_Ma.col(0) = arma::sqrt((DummyIS_u.col(0) % DummyIS_u.col(0) + DummyIS_v.col(0) % DummyIS_v.col(0)) / k / R / DummyIS_T.col(0));

	for (nd = 1; nd < n_dummy; nd++)
	{
		DummyIS_Den.col(nd) = DummyIS_Den.col(0);
		DummyIS_E.col(nd) = DummyIS_E.col(0);
		DummyIS_P.col(nd) = DummyIS_P.col(0);
		DummyIS_T.col(nd) = DummyIS_T.col(0);
		DummyIS_u.col(nd) = DummyIS_u.col(0);
		DummyIS_v.col(nd) = DummyIS_v.col(0);
		DummyIS_Ma.col(nd) = DummyIS_Ma.col(0);
	}

	//1.2远场边界的虚元
	//远场边界的虚元值由远场边界条件计算出的边界值给出
	DummyIN_Den.col(0) = FaceI_Den.col((long long)Jmax - 1);
	DummyIN_E.col(0) = FaceI_E.col((long long)Jmax - 1);
	DummyIN_P.col(0) = FaceI_P.col((long long)Jmax - 1);
	DummyIN_T.col(0) = FaceI_T.col((long long)Jmax - 1);
	DummyIN_u.col(0) = FaceI_u.col((long long)Jmax - 1);
	DummyIN_v.col(0) = FaceI_v.col((long long)Jmax - 1);
	DummyIN_Ma.col(0) = arma::sqrt((DummyIN_u.col(0) % DummyIN_u.col(0) + DummyIN_v.col(0) % DummyIN_v.col(0)) / k / R / DummyIN_T.col(0));

	for (nd = 1; nd < n_dummy; nd++)
	{
		DummyIN_Den.col(nd) = DummyIN_Den.col(0);
		DummyIN_E.col(nd) = DummyIN_E.col(0);
		DummyIN_P.col(nd) = DummyIN_P.col(0);
		DummyIN_T.col(nd) = DummyIN_T.col(0);
		DummyIN_u.col(nd) = DummyIN_u.col(0);
		DummyIN_v.col(nd) = DummyIN_v.col(0);
		DummyIN_Ma.col(nd) = DummyIN_Ma.col(0);
	}

	//1.3割缝边界的虚元
	for (nd = 0; nd < n_dummy; nd++)
	{
		DummyJW_Den.row(nd) = Cell_Den.row((long long)Imax - 2 - nd);
		DummyJW_E.row(nd) = Cell_E.row((long long)Imax - 2 - nd);
		DummyJW_P.row(nd) = Cell_P.row((long long)Imax - 2 - nd);
		DummyJW_T.row(nd) = Cell_T.row((long long)Imax - 2 - nd);
		DummyJW_u.row(nd) = Cell_u.row((long long)Imax - 2 - nd);
		DummyJW_v.row(nd) = Cell_v.row((long long)Imax - 2 - nd);
		DummyJW_Ma.row(nd) = arma::sqrt((DummyJW_u.row(nd) % DummyJW_u.row(nd) + DummyJW_v.row(nd) % DummyJW_v.row(nd)) / k / R / DummyJW_T.row(nd));
	}
	DummyJE_Den = Cell_Den.rows(0, (long long)n_dummy - 1);
	DummyJE_E = Cell_E.rows(0, (long long)n_dummy - 1);
	DummyJE_P = Cell_P.rows(0, (long long)n_dummy - 1);
	DummyJE_T = Cell_T.rows(0, (long long)n_dummy - 1);
	DummyJE_u = Cell_u.rows(0, (long long)n_dummy - 1);
	DummyJE_v = Cell_v.rows(0, (long long)n_dummy - 1);
	DummyJE_Ma = arma::sqrt((DummyJE_u % DummyJE_u + DummyJE_v % DummyJE_v) / k / R / DummyJW_T);

	//2.计算虚元守恒量
	DummyIS_U1 = DummyIS_Den;
	DummyIS_U2 = DummyIS_Den % DummyIS_u;
	DummyIS_U3 = DummyIS_Den % DummyIS_v;
	DummyIS_U4 = DummyIS_Den % DummyIS_E;

	DummyIN_U1 = DummyIN_Den;
	DummyIN_U2 = DummyIN_Den % DummyIN_u;
	DummyIN_U3 = DummyIN_Den % DummyIN_v;
	DummyIN_U4 = DummyIN_Den % DummyIN_E;

	DummyJW_U1 = DummyJW_Den;
	DummyJW_U2 = DummyJW_Den % DummyJW_u;
	DummyJW_U3 = DummyJW_Den % DummyJW_v;
	DummyJW_U4 = DummyJW_Den % DummyJW_E;

	DummyJE_U1 = DummyJE_Den;
	DummyJE_U2 = DummyJE_Den % DummyJE_u;
	DummyJE_U3 = DummyJE_Den % DummyJE_v;
	DummyJE_U4 = DummyJE_Den % DummyJE_E;
}
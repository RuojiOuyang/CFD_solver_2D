/*
DynamicArray函数用于动态分配内存
*/

#include "VariablesDeclear.h"

using namespace std;
using namespace arma;

void DynamicArray()
{
	Cell_x.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_y.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_Vol.zeros((long long)Imax - 1, (long long)Jmax - 1);
	FaceI_nx.zeros((long long)Imax - 1, Jmax);
	FaceI_ny.zeros((long long)Imax - 1, Jmax);
	FaceJ_nx.zeros(Imax, (long long)Jmax - 1);
	FaceJ_ny.zeros(Imax, (long long)Jmax - 1);
	Node_x.zeros(Imax, Jmax);
	Node_y.zeros(Imax, Jmax);

	P_former.zeros((long long)Imax - 1, (long long)Jmax - 1);
	muI.zeros((long long)Imax - 1, (long long)Jmax + 2 * (long long)n_dummy - 3);
	muJ.zeros((long long)Imax + 2 * (long long)n_dummy - 3, (long long)Jmax - 1);

	Node_Den.zeros(Imax, Jmax);
	Node_u.zeros(Imax, Jmax);
	Node_v.zeros(Imax, Jmax);
	Node_E.zeros(Imax, Jmax);
	Node_P.zeros(Imax, Jmax);
	Node_T.zeros(Imax, Jmax);
	Node_Ma.zeros(Imax, Jmax);

	//Cell是(Imax-1)*(Jmax-1)维的
	//IS是(Imax-1)*n_dummy维的
	//IN是(Imax-1)*n_dummy维的
	//JW是n_dummy*(Jmax-1)维的
	//JE是n_dummy*(Jmax-1)维的
	//FaceI是(Imax-1)*Jmax维的
	//FaceJ是Imax*(Jmax-1)维的

	FaceI_Den.zeros((long long)Imax - 1, Jmax);
	FaceI_u.zeros((long long)Imax - 1, Jmax);
	FaceI_v.zeros((long long)Imax - 1, Jmax);
	FaceI_E.zeros((long long)Imax - 1, Jmax);
	FaceI_P.zeros((long long)Imax - 1, Jmax);
	FaceI_T.zeros((long long)Imax - 1, Jmax);
	FaceI_U1.zeros((long long)Imax - 1, Jmax);
	FaceI_U2.zeros((long long)Imax - 1, Jmax);
	FaceI_U3.zeros((long long)Imax - 1, Jmax);
	FaceI_U4.zeros((long long)Imax - 1, Jmax);
	FaceI_F1.zeros((long long)Imax - 1, Jmax);
	FaceI_F2.zeros((long long)Imax - 1, Jmax);
	FaceI_F3.zeros((long long)Imax - 1, Jmax);
	FaceI_F4.zeros((long long)Imax - 1, Jmax);
	FaceI_lamta.zeros((long long)Imax - 1, Jmax);
	FaceI_epslion1.zeros((long long)Imax - 1, Jmax);
	FaceI_epslion2.zeros((long long)Imax - 1, Jmax);
	FaceI_dU1.zeros((long long)Imax - 1, Jmax);
	FaceI_dU2.zeros((long long)Imax - 1, Jmax);
	FaceI_dU3.zeros((long long)Imax - 1, Jmax);
	FaceI_dU4.zeros((long long)Imax - 1, Jmax);
	FaceI_dddU1.zeros((long long)Imax - 1, Jmax);
	FaceI_dddU2.zeros((long long)Imax - 1, Jmax);
	FaceI_dddU3.zeros((long long)Imax - 1, Jmax);
	FaceI_dddU4.zeros((long long)Imax - 1, Jmax);
	FaceI_D1.zeros((long long)Imax - 1, Jmax);
	FaceI_D2.zeros((long long)Imax - 1, Jmax);
	FaceI_D3.zeros((long long)Imax - 1, Jmax);
	FaceI_D4.zeros((long long)Imax - 1, Jmax);

	FaceJ_Den.zeros(Imax, (long long)Jmax - 1);
	FaceJ_u.zeros(Imax, (long long)Jmax - 1);
	FaceJ_v.zeros(Imax, (long long)Jmax - 1);
	FaceJ_E.zeros(Imax, (long long)Jmax - 1);
	FaceJ_P.zeros(Imax, (long long)Jmax - 1);
	FaceJ_T.zeros(Imax, (long long)Jmax - 1);
	FaceJ_U1.zeros(Imax, (long long)Jmax - 1);
	FaceJ_U2.zeros(Imax, (long long)Jmax - 1);
	FaceJ_U3.zeros(Imax, (long long)Jmax - 1);
	FaceJ_U4.zeros(Imax, (long long)Jmax - 1);
	FaceJ_F1.zeros(Imax, (long long)Jmax - 1);
	FaceJ_F2.zeros(Imax, (long long)Jmax - 1);
	FaceJ_F3.zeros(Imax, (long long)Jmax - 1);
	FaceJ_F4.zeros(Imax, (long long)Jmax - 1);
	FaceJ_lamta.zeros(Imax, (long long)Jmax - 1);
	FaceJ_epslion1.zeros(Imax, (long long)Jmax - 1);
	FaceJ_epslion2.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dU1.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dU2.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dU3.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dU4.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dddU1.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dddU2.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dddU3.zeros(Imax, (long long)Jmax - 1);
	FaceJ_dddU4.zeros(Imax, (long long)Jmax - 1);
	FaceJ_D1.zeros(Imax, (long long)Jmax - 1);
	FaceJ_D2.zeros(Imax, (long long)Jmax - 1);
	FaceJ_D3.zeros(Imax, (long long)Jmax - 1);
	FaceJ_D4.zeros(Imax, (long long)Jmax - 1);

	Cell_Den.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_u.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_v.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_E.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_P.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_T.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_Ma.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_dt.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_c.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U1.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U2.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U3.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U4.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U1_former.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U2_former.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U3_former.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_U4_former.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_QC1.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_QC2.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_QC3.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_QC4.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_D1.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_D2.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_D3.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_D4.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_Med1.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_Med2.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_Med3.zeros((long long)Imax - 1, (long long)Jmax - 1);
	Cell_Med4.zeros((long long)Imax - 1, (long long)Jmax - 1);

	DummyIS_Den.zeros((long long)Imax - 1, n_dummy);
	DummyIS_u.zeros((long long)Imax - 1, n_dummy);
	DummyIS_v.zeros((long long)Imax - 1, n_dummy);
	DummyIS_E.zeros((long long)Imax - 1, n_dummy);
	DummyIS_P.zeros((long long)Imax - 1, n_dummy);
	DummyIS_T.zeros((long long)Imax - 1, n_dummy);
	DummyIS_Ma.zeros((long long)Imax - 1, n_dummy);
	DummyIS_dt.zeros((long long)Imax - 1, n_dummy);
	DummyIS_c.zeros((long long)Imax - 1, n_dummy);
	DummyIS_U1.zeros((long long)Imax - 1, n_dummy);
	DummyIS_U2.zeros((long long)Imax - 1, n_dummy);
	DummyIS_U3.zeros((long long)Imax - 1, n_dummy);
	DummyIS_U4.zeros((long long)Imax - 1, n_dummy);

	DummyIN_Den.zeros((long long)Imax - 1, n_dummy);
	DummyIN_u.zeros((long long)Imax - 1, n_dummy);
	DummyIN_v.zeros((long long)Imax - 1, n_dummy);
	DummyIN_E.zeros((long long)Imax - 1, n_dummy);
	DummyIN_P.zeros((long long)Imax - 1, n_dummy);
	DummyIN_T.zeros((long long)Imax - 1, n_dummy);
	DummyIN_Ma.zeros((long long)Imax - 1, n_dummy);
	DummyIN_dt.zeros((long long)Imax - 1, n_dummy);
	DummyIN_c.zeros((long long)Imax - 1, n_dummy);
	DummyIN_U1.zeros((long long)Imax - 1, n_dummy);
	DummyIN_U2.zeros((long long)Imax - 1, n_dummy);
	DummyIN_U3.zeros((long long)Imax - 1, n_dummy);
	DummyIN_U4.zeros((long long)Imax - 1, n_dummy);

	DummyJW_Den.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_u.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_v.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_E.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_P.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_T.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_Ma.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_dt.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_c.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_U1.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_U2.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_U3.zeros(n_dummy, (long long)Jmax - 1);
	DummyJW_U4.zeros(n_dummy, (long long)Jmax - 1);

	DummyJE_Den.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_u.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_v.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_E.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_P.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_T.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_Ma.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_dt.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_c.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_U1.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_U2.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_U3.zeros(n_dummy, (long long)Jmax - 1);
	DummyJE_U4.zeros(n_dummy, (long long)Jmax - 1);

	//RK系数
	RKalpha = new double[6];//RK的系数
}
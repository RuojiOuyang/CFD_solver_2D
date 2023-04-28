/*
Dissipative函数计算人工黏性项，每个中心点处的Dij
Max函数计算四个数中的最大值
*/

#include "Main.h"

using namespace std;

extern double Max(double a, double b, double c, double d);

void Dissipative()
{
	//1.计算lamta
	//1.1边界上的当地通量Jacobi矩阵谱半径的某种近似
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][1].lamta = Cell[i][1].Vol / Cell[i][1].localdt * CFL;
		FaceI[i][Jmax].lamta = Cell[i][Jmax - 1].Vol / Cell[i][Jmax - 1].localdt * CFL;
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		FaceJ[1][j].lamta = 0.5 * (Cell[1][j].Vol / Cell[1][j].localdt * CFL + Cell[Imax - 1][j].Vol / Cell[Imax - 1][j].localdt * CFL);//这里这么写是因为是割缝边界
		FaceJ[Imax][j].lamta = 0.5 * (Cell[Imax - 1][j].Vol / Cell[Imax - 1][j].localdt * CFL + Cell[1][j].Vol / Cell[1][j].localdt * CFL);
	}
	//1.2内部边界上的当地通量Jacobi矩阵谱半径的某种近似
	for (j = 2; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			FaceI[i][j].lamta = 0.5 * (Cell[i][j].Vol / Cell[i][j].localdt + Cell[i][j - 1].Vol / Cell[i][j - 1].localdt) * CFL;
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 2; i <= Imax - 1; i++)
		{
			FaceJ[i][j].lamta = 0.5 * (Cell[i][j].Vol / Cell[i][j].localdt + Cell[i - 1][j].Vol / Cell[i - 1][j].localdt) * CFL;
		}
	}

	//2.激波感受因子
	//2.1边界
	for (i = 1; i <= Imax - 1; i++)
	{
		muI[i][1] = fabs(DummyIS[i][1].P - 2.0 * DummyIS[i][2].P + DummyIS[i][3].P) / fabs(DummyIS[i][1].P + 2.0 * DummyIS[i][2].P + DummyIS[i][3].P);
		muI[i][2] = fabs(Cell[i][1].P - 2.0 * DummyIS[i][1].P + DummyIS[i][2].P) / fabs(Cell[i][1].P + 2.0 * DummyIS[i][1].P + DummyIS[i][2].P);
		muI[i][3] = fabs(Cell[i][2].P - 2.0 * Cell[i][1].P + DummyIS[i][1].P) / fabs(Cell[i][2].P + 2.0 * Cell[i][1].P + DummyIS[i][1].P);
		muI[i][Jmax + 1] = fabs(DummyIN[i][1].P - 2.0 * Cell[i][Jmax - 1].P + Cell[i][Jmax - 2].P) / fabs(DummyIN[i][1].P + 2.0 * Cell[i][Jmax - 1].P + Cell[i][Jmax - 2].P);
		muI[i][Jmax + 2] = fabs(DummyIN[i][2].P - 2.0 * DummyIN[i][1].P + Cell[i][Jmax - 1].P) / fabs(DummyIN[i][2].P + 2.0 * DummyIN[i][1].P + Cell[i][Jmax - 1].P);
		muI[i][Jmax + 3] = fabs(DummyIN[i][3].P - 2.0 * DummyIN[i][2].P + DummyIN[i][1].P) / fabs(DummyIN[i][3].P + 2.0 * DummyIN[i][2].P + DummyIN[i][1].P);
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		muJ[1][j] = fabs(DummyJW[1][j].P - 2.0 * DummyJW[2][j].P + DummyJW[3][j].P) / fabs(DummyJW[1][j].P + 2.0 * DummyJW[2][j].P + DummyJW[3][j].P);
		muJ[2][j] = fabs(Cell[1][j].P - 2.0 * DummyJW[1][j].P + DummyJW[2][j].P) / fabs(Cell[1][j].P + 2.0 * DummyJW[1][j].P + DummyJW[2][j].P);
		muJ[3][j] = fabs(Cell[2][j].P - 2.0 * Cell[1][j].P + DummyJW[1][j].P) / fabs(Cell[2][j].P + 2.0 * Cell[1][j].P + DummyJW[1][j].P);
		muJ[Imax + 1][j] = fabs(DummyJE[1][j].P - 2.0 * Cell[Imax - 1][j].P + Cell[Imax - 2][j].P) / fabs(DummyJE[1][j].P + 2.0 * Cell[Imax - 1][j].P + Cell[Imax - 2][j].P);
		muJ[Imax + 2][j] = fabs(DummyJE[2][j].P - 2.0 * DummyJE[1][j].P + Cell[Imax - 1][j].P) / fabs(DummyJE[2][j].P + 2.0 * DummyJE[1][j].P + Cell[Imax - 1][j].P);
		muJ[Imax + 3][j] = fabs(DummyJE[3][j].P - 2.0 * DummyJE[2][j].P + DummyJE[1][j].P) / fabs(DummyJE[3][j].P + 2.0 * DummyJE[2][j].P + DummyJE[1][j].P);
	}
	//2.2内部
	/*
	 _ _ _ _ _
	|_|_|_|_|_|
	|*|*|*|*|*|
	|*|*|*|*|*|
	|*|*|*|*|*|
	|_|_|_|_|_|

	*/
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 4; j <= Jmax; j++)
		{
			muI[i][j] = fabs(Cell[i][j - 1].P - 2.0 * Cell[i][j - 2].P + Cell[i][j - 3].P) / fabs(Cell[i][j - 1].P + 2.0 * Cell[i][j - 2].P + Cell[i][j - 3].P);
		}
	}
	/*
	 _ _ _ _ _
	|_|*|*|*|_|
	|_|*|*|*|_|
	|_|*|*|*|_|
	|_|*|*|*|_|
	|_|*|*|*|_|

	*/
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 4; i <= Imax; i++)
		{
			muJ[i][j] = fabs(Cell[i - 1][j].P - 2.0 * Cell[i - 2][j].P + Cell[i - 3][j].P) / fabs(Cell[i - 1][j].P + 2.0 * Cell[i - 2][j].P + Cell[i - 3][j].P);
		}
	}

	//3.自适应粘性系数ε2和ε4
	double k2, k4;//ε的比例系数
	k2 = 0.5;
	k4 = 1.0 / 128.0;
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax; j++)
		{
			FaceI[i][j].epslion[1] = k2 * Max(muI[i][j + 3], muI[i][j + 2], muI[i][j + 1], muI[i][j]);
			FaceI[i][j].epslion[2] = max(0.0, k4 - FaceI[i][j].epslion[1]);
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax; i++)
		{
			FaceJ[i][j].epslion[1] = k2 * Max(muJ[i + 3][j], muJ[i + 2][j], muJ[i + 1][j], muJ[i][j]);
			FaceJ[i][j].epslion[2] = max(0.0, k4 - FaceJ[i][j].epslion[1]);
		}
	}

	//4.计算差分
	//4.1边界
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][1].dU[1] = Cell[i][1].U[1] - DummyIS[i][1].U[1];
		FaceI[i][1].dU[2] = Cell[i][1].U[2] - DummyIS[i][1].U[2];
		FaceI[i][1].dU[3] = Cell[i][1].U[3] - DummyIS[i][1].U[3];
		FaceI[i][1].dU[4] = Cell[i][1].U[4] - DummyIS[i][1].U[4];
		FaceI[i][Jmax].dU[1] = DummyIN[i][1].U[1] - Cell[i][Jmax - 1].U[1];
		FaceI[i][Jmax].dU[2] = DummyIN[i][1].U[2] - Cell[i][Jmax - 1].U[2];
		FaceI[i][Jmax].dU[3] = DummyIN[i][1].U[3] - Cell[i][Jmax - 1].U[3];
		FaceI[i][Jmax].dU[4] = DummyIN[i][1].U[4] - Cell[i][Jmax - 1].U[4];
		FaceI[i][1].dddU[1] = Cell[i][2].U[1] - 3.0 * Cell[i][1].U[1] + 3.0 * DummyIS[i][1].U[1] - DummyIS[i][2].U[1];
		FaceI[i][1].dddU[2] = Cell[i][2].U[2] - 3.0 * Cell[i][1].U[2] + 3.0 * DummyIS[i][1].U[2] - DummyIS[i][2].U[2];
		FaceI[i][1].dddU[3] = Cell[i][2].U[3] - 3.0 * Cell[i][1].U[3] + 3.0 * DummyIS[i][1].U[3] - DummyIS[i][2].U[3];
		FaceI[i][1].dddU[4] = Cell[i][2].U[4] - 3.0 * Cell[i][1].U[4] + 3.0 * DummyIS[i][1].U[4] - DummyIS[i][2].U[4];
		FaceI[i][2].dddU[1] = Cell[i][3].U[1] - 3.0 * Cell[i][2].U[1] + 3.0 * Cell[i][1].U[1] - DummyIS[i][1].U[1];
		FaceI[i][2].dddU[2] = Cell[i][3].U[2] - 3.0 * Cell[i][2].U[2] + 3.0 * Cell[i][1].U[2] - DummyIS[i][1].U[2];
		FaceI[i][2].dddU[3] = Cell[i][3].U[3] - 3.0 * Cell[i][2].U[3] + 3.0 * Cell[i][1].U[3] - DummyIS[i][1].U[3];
		FaceI[i][2].dddU[4] = Cell[i][3].U[4] - 3.0 * Cell[i][2].U[4] + 3.0 * Cell[i][1].U[4] - DummyIS[i][1].U[4];
		FaceI[i][Jmax].dddU[1] = DummyIN[i][2].U[1] - 3.0 * DummyIN[i][1].U[1] + 3.0 * Cell[i][Jmax - 1].U[1] - Cell[i][Jmax - 2].U[1];
		FaceI[i][Jmax].dddU[2] = DummyIN[i][2].U[2] - 3.0 * DummyIN[i][1].U[2] + 3.0 * Cell[i][Jmax - 1].U[2] - Cell[i][Jmax - 2].U[2];
		FaceI[i][Jmax].dddU[3] = DummyIN[i][2].U[3] - 3.0 * DummyIN[i][1].U[3] + 3.0 * Cell[i][Jmax - 1].U[3] - Cell[i][Jmax - 2].U[3];
		FaceI[i][Jmax].dddU[4] = DummyIN[i][2].U[4] - 3.0 * DummyIN[i][1].U[4] + 3.0 * Cell[i][Jmax - 1].U[4] - Cell[i][Jmax - 2].U[4];
		FaceI[i][Jmax - 1].dddU[1] = DummyIN[i][1].U[1] - 3.0 * Cell[i][Jmax - 1].U[1] + 3.0 * Cell[i][Jmax - 2].U[1] - Cell[i][Jmax - 3].U[1];
		FaceI[i][Jmax - 1].dddU[2] = DummyIN[i][1].U[2] - 3.0 * Cell[i][Jmax - 1].U[2] + 3.0 * Cell[i][Jmax - 2].U[2] - Cell[i][Jmax - 3].U[2];
		FaceI[i][Jmax - 1].dddU[3] = DummyIN[i][1].U[3] - 3.0 * Cell[i][Jmax - 1].U[3] + 3.0 * Cell[i][Jmax - 2].U[3] - Cell[i][Jmax - 3].U[3];
		FaceI[i][Jmax - 1].dddU[4] = DummyIN[i][1].U[4] - 3.0 * Cell[i][Jmax - 1].U[4] + 3.0 * Cell[i][Jmax - 2].U[4] - Cell[i][Jmax - 3].U[4];
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		FaceJ[1][j].dU[1] = Cell[1][j].U[1] - DummyJW[1][j].U[1];
		FaceJ[1][j].dU[2] = Cell[1][j].U[2] - DummyJW[1][j].U[2];
		FaceJ[1][j].dU[3] = Cell[1][j].U[3] - DummyJW[1][j].U[3];
		FaceJ[1][j].dU[4] = Cell[1][j].U[4] - DummyJW[1][j].U[4];
		FaceJ[Imax][j].dU[1] = DummyJE[1][j].U[1] - Cell[Imax - 1][j].U[1];
		FaceJ[Imax][j].dU[2] = DummyJE[1][j].U[2] - Cell[Imax - 1][j].U[2];
		FaceJ[Imax][j].dU[3] = DummyJE[1][j].U[3] - Cell[Imax - 1][j].U[3];
		FaceJ[Imax][j].dU[4] = DummyJE[1][j].U[4] - Cell[Imax - 1][j].U[4];
		FaceJ[1][j].dddU[1] = Cell[2][j].U[1] - 3.0 * Cell[1][j].U[1] + 3.0 * DummyJW[1][j].U[1] - DummyJW[2][j].U[1];
		FaceJ[1][j].dddU[2] = Cell[2][j].U[2] - 3.0 * Cell[1][j].U[2] + 3.0 * DummyJW[1][j].U[2] - DummyJW[2][j].U[2];
		FaceJ[1][j].dddU[3] = Cell[2][j].U[3] - 3.0 * Cell[1][j].U[3] + 3.0 * DummyJW[1][j].U[3] - DummyJW[2][j].U[3];
		FaceJ[1][j].dddU[4] = Cell[2][j].U[4] - 3.0 * Cell[1][j].U[4] + 3.0 * DummyJW[1][j].U[4] - DummyJW[2][j].U[4];
		FaceJ[2][j].dddU[1] = Cell[3][j].U[1] - 3.0 * Cell[2][j].U[1] + 3.0 * Cell[1][j].U[1] - DummyJW[1][j].U[1];
		FaceJ[2][j].dddU[2] = Cell[3][j].U[2] - 3.0 * Cell[2][j].U[2] + 3.0 * Cell[1][j].U[2] - DummyJW[1][j].U[2];
		FaceJ[2][j].dddU[3] = Cell[3][j].U[3] - 3.0 * Cell[2][j].U[3] + 3.0 * Cell[1][j].U[3] - DummyJW[1][j].U[3];
		FaceJ[2][j].dddU[4] = Cell[3][j].U[4] - 3.0 * Cell[2][j].U[4] + 3.0 * Cell[1][j].U[4] - DummyJW[1][j].U[4];
		FaceJ[Imax][j].dddU[1] = DummyJE[2][j].U[1] - 3.0 * DummyJE[1][j].U[1] + 3.0 * Cell[Imax - 1][j].U[1] - Cell[Imax - 2][j].U[1];
		FaceJ[Imax][j].dddU[2] = DummyJE[2][j].U[2] - 3.0 * DummyJE[1][j].U[2] + 3.0 * Cell[Imax - 1][j].U[2] - Cell[Imax - 2][j].U[2];
		FaceJ[Imax][j].dddU[3] = DummyJE[2][j].U[3] - 3.0 * DummyJE[1][j].U[3] + 3.0 * Cell[Imax - 1][j].U[3] - Cell[Imax - 2][j].U[3];
		FaceJ[Imax][j].dddU[4] = DummyJE[2][j].U[4] - 3.0 * DummyJE[1][j].U[4] + 3.0 * Cell[Imax - 1][j].U[4] - Cell[Imax - 2][j].U[4];
		FaceJ[Imax - 1][j].dddU[1] = DummyJE[1][j].U[1] - 3.0 * Cell[Imax - 1][j].U[1] + 3.0 * Cell[Imax - 2][j].U[1] - Cell[Imax - 3][j].U[1];
		FaceJ[Imax - 1][j].dddU[2] = DummyJE[1][j].U[2] - 3.0 * Cell[Imax - 1][j].U[2] + 3.0 * Cell[Imax - 2][j].U[2] - Cell[Imax - 3][j].U[2];
		FaceJ[Imax - 1][j].dddU[3] = DummyJE[1][j].U[3] - 3.0 * Cell[Imax - 1][j].U[3] + 3.0 * Cell[Imax - 2][j].U[3] - Cell[Imax - 3][j].U[3];
		FaceJ[Imax - 1][j].dddU[4] = DummyJE[1][j].U[4] - 3.0 * Cell[Imax - 1][j].U[4] + 3.0 * Cell[Imax - 2][j].U[4] - Cell[Imax - 3][j].U[4];
	}

	//4.2内部
	//4.2.1内部一阶差分
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 2; j <= Jmax - 1; j++)
		{
			FaceI[i][j].dU[1] = Cell[i][j].U[1] - Cell[i][j - 1].U[1];
			FaceI[i][j].dU[2] = Cell[i][j].U[2] - Cell[i][j - 1].U[2];
			FaceI[i][j].dU[3] = Cell[i][j].U[3] - Cell[i][j - 1].U[3];
			FaceI[i][j].dU[4] = Cell[i][j].U[4] - Cell[i][j - 1].U[4];
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 2; i <= Imax - 1; i++)
		{
			FaceJ[i][j].dU[1] = Cell[i][j].U[1] - Cell[i - 1][j].U[1];
			FaceJ[i][j].dU[2] = Cell[i][j].U[2] - Cell[i - 1][j].U[2];
			FaceJ[i][j].dU[3] = Cell[i][j].U[3] - Cell[i - 1][j].U[3];
			FaceJ[i][j].dU[4] = Cell[i][j].U[4] - Cell[i - 1][j].U[4];
		}
	}
	//4.2.2内部三阶差分
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 3; j <= Jmax - 2; j++)
		{
			FaceI[i][j].dddU[1] = Cell[i][j + 1].U[1] - 3.0 * Cell[i][j].U[1] + 3.0 * Cell[i][j - 1].U[1] - Cell[i][j - 2].U[1];
			FaceI[i][j].dddU[2] = Cell[i][j + 1].U[2] - 3.0 * Cell[i][j].U[2] + 3.0 * Cell[i][j - 1].U[2] - Cell[i][j - 2].U[2];
			FaceI[i][j].dddU[3] = Cell[i][j + 1].U[3] - 3.0 * Cell[i][j].U[3] + 3.0 * Cell[i][j - 1].U[3] - Cell[i][j - 2].U[3];
			FaceI[i][j].dddU[4] = Cell[i][j + 1].U[4] - 3.0 * Cell[i][j].U[4] + 3.0 * Cell[i][j - 1].U[4] - Cell[i][j - 2].U[4];
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 3; i <= Imax - 2; i++)
		{
			FaceJ[i][j].dddU[1] = Cell[i + 1][j].U[1] - 3.0 * Cell[i][j].U[1] + 3.0 * Cell[i - 1][j].U[1] - Cell[i - 2][j].U[1];
			FaceJ[i][j].dddU[2] = Cell[i + 1][j].U[2] - 3.0 * Cell[i][j].U[2] + 3.0 * Cell[i - 1][j].U[2] - Cell[i - 2][j].U[2];
			FaceJ[i][j].dddU[3] = Cell[i + 1][j].U[3] - 3.0 * Cell[i][j].U[3] + 3.0 * Cell[i - 1][j].U[3] - Cell[i - 2][j].U[3];
			FaceJ[i][j].dddU[4] = Cell[i + 1][j].U[4] - 3.0 * Cell[i][j].U[4] + 3.0 * Cell[i - 1][j].U[4] - Cell[i - 2][j].U[4];
		}
	}

	//5.所有边界上的人工黏性项
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax; j++)
		{
			FaceI[i][j].FD[1] = FaceI[i][j].lamta * (FaceI[i][j].epslion[1] * FaceI[i][j].dU[1] - FaceI[i][j].epslion[2] * FaceI[i][j].dddU[1]);
			FaceI[i][j].FD[2] = FaceI[i][j].lamta * (FaceI[i][j].epslion[1] * FaceI[i][j].dU[2] - FaceI[i][j].epslion[2] * FaceI[i][j].dddU[2]);
			FaceI[i][j].FD[3] = FaceI[i][j].lamta * (FaceI[i][j].epslion[1] * FaceI[i][j].dU[3] - FaceI[i][j].epslion[2] * FaceI[i][j].dddU[3]);
			FaceI[i][j].FD[4] = FaceI[i][j].lamta * (FaceI[i][j].epslion[1] * FaceI[i][j].dU[4] - FaceI[i][j].epslion[2] * FaceI[i][j].dddU[4]);
		}
	}
	for (i = 1; i <= Imax; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			FaceJ[i][j].FD[1] = FaceJ[i][j].lamta * (FaceJ[i][j].epslion[1] * FaceJ[i][j].dU[1] - FaceJ[i][j].epslion[2] * FaceJ[i][j].dddU[1]);
			FaceJ[i][j].FD[2] = FaceJ[i][j].lamta * (FaceJ[i][j].epslion[1] * FaceJ[i][j].dU[2] - FaceJ[i][j].epslion[2] * FaceJ[i][j].dddU[2]);
			FaceJ[i][j].FD[3] = FaceJ[i][j].lamta * (FaceJ[i][j].epslion[1] * FaceJ[i][j].dU[3] - FaceJ[i][j].epslion[2] * FaceJ[i][j].dddU[3]);
			FaceJ[i][j].FD[4] = FaceJ[i][j].lamta * (FaceJ[i][j].epslion[1] * FaceJ[i][j].dU[4] - FaceJ[i][j].epslion[2] * FaceJ[i][j].dddU[4]);
		}
	}

	//6.控制体人工黏性
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			Cell[i][j].D[1] = FaceI[i][j + 1].FD[1] - FaceI[i][j].FD[1] + FaceJ[i + 1][j].FD[1] - FaceJ[i][j].FD[1];
			Cell[i][j].D[2] = FaceI[i][j + 1].FD[2] - FaceI[i][j].FD[2] + FaceJ[i + 1][j].FD[2] - FaceJ[i][j].FD[2];
			Cell[i][j].D[3] = FaceI[i][j + 1].FD[3] - FaceI[i][j].FD[3] + FaceJ[i + 1][j].FD[3] - FaceJ[i][j].FD[3];
			Cell[i][j].D[4] = FaceI[i][j + 1].FD[4] - FaceI[i][j].FD[4] + FaceJ[i + 1][j].FD[4] - FaceJ[i][j].FD[4];
		}
	}
}

double Max(double a, double b, double c, double d)
{
	double e;
	e = a;
	if (a <= b)
	{
		e = b;
	}
	if (e <= c)
	{
		e = c;
	}
	if (e <= d)
	{
		e = d;
	}
	return e;
}
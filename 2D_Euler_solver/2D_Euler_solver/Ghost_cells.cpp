/*
Ghost_cells函数，处理虚元，为人工黏性项的计算做铺垫
*/

#include "Main.h"

using namespace std;

void Ghost_cells()
{
	int nd;//本函数中需要用到的迭代标记，局部变量

	//1.算出各个边界处的虚元
	//1.1物面边界的虚元
	double V_f, V_a, V_b;
	for (i = 1; i <= Imax - 1; i++)
	{
		DummyIS[i][1].Den = Cell[i][1].Den;
		DummyIS[i][1].E = Cell[i][1].E;
		DummyIS[i][1].P = Cell[i][1].P;
		DummyIS[i][1].T = Cell[i][1].T;
		DummyIS[i][1].H = Cell[i][1].H;
		V_f = Cell[i][1].u * FaceI[i][1].nx + Cell[i][1].v * FaceI[i][1].ny;
		V_a = V_f * FaceI[i][1].nx / (FaceI[i][1].nx * FaceI[i][1].nx + FaceI[i][1].ny * FaceI[i][1].ny);
		V_b = V_f * FaceI[i][1].ny / (FaceI[i][1].nx * FaceI[i][1].nx + FaceI[i][1].ny * FaceI[i][1].ny);
		DummyIS[i][1].u = Cell[i][1].u - 2.0 * V_a;
		DummyIS[i][1].v = Cell[i][1].v - 2.0 * V_b;
		DummyIS[i][1].Ma = sqrt((DummyIS[i][1].u * DummyIS[i][1].u + DummyIS[i][1].v * DummyIS[i][1].v) / k / R / DummyIS[i][1].T);

		for (nd = 2; nd <= n_dummy; nd++)
		{
			DummyIS[i][nd] = DummyIS[i][1];
		}
	}

	//1.2远场边界的虚元
	//远场边界的虚元值由远场边界条件计算出的边界值给出
	for (i = 1; i <= Imax - 1; i++)
	{
		DummyIN[i][1].Den = FaceI[i][Jmax].Den;
		DummyIN[i][1].E = FaceI[i][Jmax].E;
		DummyIN[i][1].P = FaceI[i][Jmax].P;
		DummyIN[i][1].T = FaceI[i][Jmax].T;
		DummyIN[i][1].H = FaceI[i][Jmax].E + FaceI[i][Jmax].P / FaceI[i][Jmax].Den;
		DummyIN[i][1].u = FaceI[i][Jmax].u;
		DummyIN[i][1].v = FaceI[i][Jmax].v;
		DummyIN[i][1].Ma = sqrt((DummyIN[i][1].u * DummyIN[i][1].u + DummyIN[i][1].v * DummyIN[i][1].v) / k / R / DummyIN[i][1].T);

		for (nd = 2; nd <= n_dummy; nd++)
		{
			DummyIN[i][nd] = DummyIN[i][1];
		}
	}

	//1.3割缝边界的虚元
	for (j = 1; j <= Jmax - 1; j++)//左边界
	{
		for ( nd = 1; nd <= n_dummy; nd++)
		{
			//cout << DummyJW[nd][j].Den << endl;
			DummyJW[nd][j].Den = Cell[Imax - nd][j].Den;
			DummyJW[nd][j].E = Cell[Imax - nd][j].E;
			DummyJW[nd][j].P = Cell[Imax - nd][j].P;
			DummyJW[nd][j].T = Cell[Imax - nd][j].T;
			DummyJW[nd][j].H = Cell[Imax - nd][j].H;
			DummyJW[nd][j].u = Cell[Imax - nd][j].u;
			DummyJW[nd][j].v = Cell[Imax - nd][j].v;
			DummyJW[nd][j].Ma = sqrt((DummyJW[nd][j].u * DummyJW[nd][j].u + DummyJW[nd][j].v * DummyJW[nd][j].v) / k / R / DummyJW[nd][j].T);
		}
	}
	for (j = 1; j <= Jmax - 1; j++)//右边界
	{
		for (nd = 1; nd <= n_dummy; nd++)
		{
			DummyJE[nd][j].Den = Cell[nd][j].Den;
			DummyJE[nd][j].E = Cell[nd][j].E;
			DummyJE[nd][j].P = Cell[nd][j].P;
			DummyJE[nd][j].T = Cell[nd][j].T;
			DummyJE[nd][j].H = Cell[nd][j].H;
			DummyJE[nd][j].u = Cell[nd][j].u;
			DummyJE[nd][j].v = Cell[nd][j].v;
			DummyJE[nd][j].Ma = sqrt((DummyJE[nd][j].u * DummyJE[nd][j].u + DummyJE[nd][j].v * DummyJE[nd][j].v) / k / R / DummyJE[nd][j].T);
		}
	}

	//2.计算虚元守恒量
	for (nd = 1; nd <= n_dummy; nd++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			DummyIS[i][nd].U[1] = DummyIS[i][nd].Den;
			DummyIS[i][nd].U[2] = DummyIS[i][nd].Den * DummyIS[i][nd].u;
			DummyIS[i][nd].U[3] = DummyIS[i][nd].Den * DummyIS[i][nd].v;
			DummyIS[i][nd].U[4] = DummyIS[i][nd].Den * DummyIS[i][nd].E;

			DummyIN[i][nd].U[1] = DummyIN[i][nd].Den;
			DummyIN[i][nd].U[2] = DummyIN[i][nd].Den * DummyIN[i][nd].u;
			DummyIN[i][nd].U[3] = DummyIN[i][nd].Den * DummyIN[i][nd].v;
			DummyIN[i][nd].U[4] = DummyIN[i][nd].Den * DummyIN[i][nd].E;
		}
		for (j = 1; j <= Jmax - 1; j++)
		{
			DummyJW[nd][j].U[1] = DummyJW[nd][j].Den;
			DummyJW[nd][j].U[2] = DummyJW[nd][j].Den * DummyJW[nd][j].u;
			DummyJW[nd][j].U[3] = DummyJW[nd][j].Den * DummyJW[nd][j].v;
			DummyJW[nd][j].U[4] = DummyJW[nd][j].Den * DummyJW[nd][j].E;

			DummyJE[nd][j].U[1] = DummyJE[nd][j].Den;
			DummyJE[nd][j].U[2] = DummyJE[nd][j].Den * DummyJE[nd][j].u;
			DummyJE[nd][j].U[3] = DummyJE[nd][j].Den * DummyJE[nd][j].v;
			DummyJE[nd][j].U[4] = DummyJE[nd][j].Den * DummyJE[nd][j].E;
		}
	}
}
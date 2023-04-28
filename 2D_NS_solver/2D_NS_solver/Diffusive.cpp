/*
Diffusive函数计算扩散项
*/

#include "Main.h"
#include "Calculate_P.h"

using namespace std;

void Diffusive()
{
	//1.计算出各个点上的Touxx等
	//1.1T、u、v在xy方向上的偏导数
	Calculate_P();
	//1.2计算Tou
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			Cell[i][j].miu = mu0 * pow(Cell[i][j].T / T0, 1.5) * (T0 + Ts) / (Cell[i][j].T + Ts);//当地气体动力黏度

			Cell[i][j].Touxx = Cell[i][j].miu * (4.0 / 3.0 * Cell[i][j].u_x - 2.0 / 3.0 * Cell[i][j].v_y);
			Cell[i][j].Touxy = Cell[i][j].miu * (Cell[i][j].u_y + Cell[i][j].v_x);
			Cell[i][j].Touyy = Cell[i][j].miu * (4.0 / 3.0 * Cell[i][j].v_y - 2.0 / 3.0 * Cell[i][j].u_x);
			Cell[i][j].Touhx = Cell[i][j].u * Cell[i][j].Touxx + Cell[i][j].v * Cell[i][j].Touxy - (-Cell[i][j].miu * cp / Pr * Cell[i][j].T_x);
			Cell[i][j].Touhy = Cell[i][j].u * Cell[i][j].Touxy + Cell[i][j].v * Cell[i][j].Touyy - (-Cell[i][j].miu * cp / Pr * Cell[i][j].T_y);
		}
	}
	//1.3第一层虚元上的Tou
	for (i = 1; i <= Imax - 1; i++)
	{
		DummyIS[i][1].miu = mu0 * pow(DummyIS[i][1].T / T0, 1.5) * (T0 + Ts) / (DummyIS[i][1].T + Ts);//当地气体动力黏度
		DummyIS[i][1].Touxx = DummyIS[i][1].miu * (4.0 / 3.0 * DummyIS[i][1].u_x - 2.0 / 3.0 * DummyIS[i][1].v_y);
		DummyIS[i][1].Touxy = DummyIS[i][1].miu * (DummyIS[i][1].u_y + DummyIS[i][1].v_x);
		DummyIS[i][1].Touyy = DummyIS[i][1].miu * (4.0 / 3.0 * DummyIS[i][1].v_y - 2.0 / 3.0 * DummyIS[i][1].u_x);
		DummyIS[i][1].Touhx = DummyIS[i][1].u * DummyIS[i][1].Touxx + DummyIS[i][1].v * DummyIS[i][1].Touxy - (-DummyIS[i][1].miu * cp / Pr * DummyIS[i][1].T_x);
		DummyIS[i][1].Touhy = DummyIS[i][1].u * DummyIS[i][1].Touxy + DummyIS[i][1].v * DummyIS[i][1].Touyy - (-DummyIS[i][1].miu * cp / Pr * DummyIS[i][1].T_y);

		DummyIN[i][1].miu = mu0 * pow(DummyIN[i][1].T / T0, 1.5) * (T0 + Ts) / (DummyIN[i][1].T + Ts);//当地气体动力黏度
		DummyIN[i][1].Touxx = DummyIN[i][1].miu * (4.0 / 3.0 * DummyIN[i][1].u_x - 2.0 / 3.0 * DummyIN[i][1].v_y);
		DummyIN[i][1].Touxy = DummyIN[i][1].miu * (DummyIN[i][1].u_y + DummyIN[i][1].v_x);
		DummyIN[i][1].Touyy = DummyIN[i][1].miu * (4.0 / 3.0 * DummyIN[i][1].v_y - 2.0 / 3.0 * DummyIN[i][1].u_x);
		DummyIN[i][1].Touhx = DummyIN[i][1].u * DummyIN[i][1].Touxx + DummyIN[i][1].v * DummyIN[i][1].Touxy - (-DummyIN[i][1].miu * cp / Pr * DummyIN[i][1].T_x);
		DummyIN[i][1].Touhy = DummyIN[i][1].u * DummyIN[i][1].Touxy + DummyIN[i][1].v * DummyIN[i][1].Touyy - (-DummyIN[i][1].miu * cp / Pr * DummyIN[i][1].T_y);
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		DummyJW[1][j].miu = mu0 * pow(DummyJW[1][j].T / T0, 1.5) * (T0 + Ts) / (DummyJW[1][j].T + Ts);//当地气体动力黏度
		DummyJW[1][j].Touxx = DummyJW[1][j].miu * (4.0 / 3.0 * DummyJW[1][j].u_x - 2.0 / 3.0 * DummyJW[1][j].v_y);
		DummyJW[1][j].Touxy = DummyJW[1][j].miu * (DummyJW[1][j].u_y + DummyJW[1][j].v_x);
		DummyJW[1][j].Touyy = DummyJW[1][j].miu * (4.0 / 3.0 * DummyJW[1][j].v_y - 2.0 / 3.0 * DummyJW[1][j].u_x);
		DummyJW[1][j].Touhx = DummyJW[1][j].u * DummyJW[1][j].Touxx + DummyJW[1][j].v * DummyJW[1][j].Touxy - (-DummyJW[1][j].miu * cp / Pr * DummyJW[1][j].T_x);
		DummyJW[1][j].Touhy = DummyJW[1][j].u * DummyJW[1][j].Touxy + DummyJW[1][j].v * DummyJW[1][j].Touyy - (-DummyJW[1][j].miu * cp / Pr * DummyJW[1][j].T_y);

		DummyJE[1][j].miu = mu0 * pow(DummyJE[1][j].T / T0, 1.5) * (T0 + Ts) / (DummyJE[1][j].T + Ts);//当地气体动力黏度
		DummyJE[1][j].Touxx = DummyJE[1][j].miu * (4.0 / 3.0 * DummyJE[1][j].u_x - 2.0 / 3.0 * DummyJE[1][j].v_y);
		DummyJE[1][j].Touxy = DummyJE[1][j].miu * (DummyJE[1][j].u_y + DummyJE[1][j].v_x);
		DummyJE[1][j].Touyy = DummyJE[1][j].miu * (4.0 / 3.0 * DummyJE[1][j].v_y - 2.0 / 3.0 * DummyJE[1][j].u_x);
		DummyJE[1][j].Touhx = DummyJE[1][j].u * DummyJE[1][j].Touxx + DummyJE[1][j].v * DummyJE[1][j].Touxy - (-DummyJE[1][j].miu * cp / Pr * DummyJE[1][j].T_x);
		DummyJE[1][j].Touhy = DummyJE[1][j].u * DummyJE[1][j].Touxy + DummyJE[1][j].v * DummyJE[1][j].Touyy - (-DummyJE[1][j].miu * cp / Pr * DummyJE[1][j].T_y);
	}

	//2.计算出所有边界上的FK
	//2.1外部边界
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][1].FK[1] = 0;
		FaceI[i][1].FK[2] = 0.5 * (Cell[i][1].Touxx + DummyIS[i][1].Touxx) * FaceI[i][1].nx + 0.5 * (Cell[i][1].Touxy + DummyIS[i][1].Touxy) * FaceI[i][1].ny;
		FaceI[i][1].FK[3] = 0.5 * (Cell[i][1].Touxy + DummyIS[i][1].Touxy) * FaceI[i][1].nx + 0.5 * (Cell[i][1].Touyy + DummyIS[i][1].Touyy) * FaceI[i][1].ny;
		FaceI[i][1].FK[4] = 0.5 * (Cell[i][1].Touhx + DummyIS[i][1].Touhx) * FaceI[i][1].nx + 0.5 * (Cell[i][1].Touhy + DummyIS[i][1].Touhy) * FaceI[i][1].ny;

		FaceI[i][Jmax].FK[1] = 0;
		FaceI[i][Jmax].FK[2] = 0.5 * (Cell[i][Jmax - 1].Touxx + DummyIN[i][1].Touxx) * FaceI[i][Jmax].nx + 0.5 * (Cell[i][Jmax - 1].Touxy + DummyIN[i][1].Touxy) * FaceI[i][Jmax].ny;
		FaceI[i][Jmax].FK[3] = 0.5 * (Cell[i][Jmax - 1].Touxy + DummyIN[i][1].Touxy) * FaceI[i][Jmax].nx + 0.5 * (Cell[i][Jmax - 1].Touyy + DummyIN[i][1].Touyy) * FaceI[i][Jmax].ny;
		FaceI[i][Jmax].FK[4] = 0.5 * (Cell[i][Jmax - 1].Touhx + DummyIN[i][1].Touhx) * FaceI[i][Jmax].nx + 0.5 * (Cell[i][Jmax - 1].Touhy + DummyIN[i][1].Touhy) * FaceI[i][Jmax].ny;
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		FaceJ[1][j].FK[1] = 0;
		FaceJ[1][j].FK[2] = 0.5 * (Cell[1][j].Touxx + DummyJW[1][j].Touxx) * FaceJ[1][j].nx + 0.5 * (Cell[1][j].Touxy + DummyJW[1][j].Touxy) * FaceJ[1][j].ny;
		FaceJ[1][j].FK[3] = 0.5 * (Cell[1][j].Touxy + DummyJW[1][j].Touxy) * FaceJ[1][j].nx + 0.5 * (Cell[1][j].Touyy + DummyJW[1][j].Touyy) * FaceJ[1][j].ny;
		FaceJ[1][j].FK[4] = 0.5 * (Cell[1][j].Touhx + DummyJW[1][j].Touhx) * FaceJ[1][j].nx + 0.5 * (Cell[1][j].Touhy + DummyJW[1][j].Touhy) * FaceJ[1][j].ny;

		FaceJ[Imax][j].FK[1] = 0;
		FaceJ[Imax][j].FK[2] = 0.5 * (DummyJE[1][j].Touxx + Cell[Imax - 1][j].Touxx) * FaceJ[Imax][j].nx + 0.5 * (DummyJE[1][j].Touxy + Cell[Imax - 1][j].Touxy) * FaceJ[Imax][j].ny;
		FaceJ[Imax][j].FK[3] = 0.5 * (DummyJE[1][j].Touxy + Cell[Imax - 1][j].Touxy) * FaceJ[Imax][j].nx + 0.5 * (DummyJE[1][j].Touyy + Cell[Imax - 1][j].Touyy) * FaceJ[Imax][j].ny;
		FaceJ[Imax][j].FK[4] = 0.5 * (DummyJE[1][j].Touhx + Cell[Imax - 1][j].Touhx) * FaceJ[Imax][j].nx + 0.5 * (DummyJE[1][j].Touhy + Cell[Imax - 1][j].Touhy) * FaceJ[Imax][j].ny;
	}
	//2.2内部边界
	for (i = 1 ; i <= Imax - 1; i++)
	{
		for (j = 2; j <= Jmax - 1; j++)
		{
			FaceI[i][j].FK[1] = 0;
			FaceI[i][j].FK[2] = 0.5 * (Cell[i][j].Touxx + Cell[i][j - 1].Touxx) * FaceI[i][j].nx + 0.5 * (Cell[i][j].Touxy + Cell[i][j - 1].Touxy) * FaceI[i][j].ny;
			FaceI[i][j].FK[3] = 0.5 * (Cell[i][j].Touxy + Cell[i][j - 1].Touxy) * FaceI[i][j].nx + 0.5 * (Cell[i][j].Touyy + Cell[i][j - 1].Touyy) * FaceI[i][j].ny;
			FaceI[i][j].FK[4] = 0.5 * (Cell[i][j].Touhx + Cell[i][j - 1].Touhx) * FaceI[i][j].nx + 0.5 * (Cell[i][j].Touhy + Cell[i][j - 1].Touhy) * FaceI[i][j].ny;
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 2; i <= Imax - 1; i++)
		{
			FaceJ[i][j].FK[1] = 0;
			FaceJ[i][j].FK[2] = 0.5 * (Cell[i][j].Touxx + Cell[i - 1][j].Touxx) * FaceJ[i][j].nx + 0.5 * (Cell[i][j].Touxy + Cell[i - 1][j].Touxy) * FaceJ[i][j].ny;
			FaceJ[i][j].FK[3] = 0.5 * (Cell[i][j].Touxy + Cell[i - 1][j].Touxy) * FaceJ[i][j].nx + 0.5 * (Cell[i][j].Touyy + Cell[i - 1][j].Touyy) * FaceJ[i][j].ny;
			FaceJ[i][j].FK[4] = 0.5 * (Cell[i][j].Touhx + Cell[i - 1][j].Touhx) * FaceJ[i][j].nx + 0.5 * (Cell[i][j].Touhy + Cell[i - 1][j].Touhy) * FaceJ[i][j].ny;
			//cout << FaceJ[i][j].FK[1] << " " << FaceJ[i][j].FK[2] << " " << FaceJ[i][j].FK[3] << " " << FaceJ[i][j].FK[4] << endl;
		}
	}

	//3.有了所有边界的FK后，控制体单元的K，也就是Qij
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].QV[1] = FaceI[i][j + 1].FK[1] - FaceI[i][j].FK[1] + FaceJ[i + 1][j].FK[1] - FaceJ[i][j].FK[1];
			Cell[i][j].QV[2] = FaceI[i][j + 1].FK[2] - FaceI[i][j].FK[2] + FaceJ[i + 1][j].FK[2] - FaceJ[i][j].FK[2];
			Cell[i][j].QV[3] = FaceI[i][j + 1].FK[3] - FaceI[i][j].FK[3] + FaceJ[i + 1][j].FK[3] - FaceJ[i][j].FK[3];
			Cell[i][j].QV[4] = FaceI[i][j + 1].FK[4] - FaceI[i][j].FK[4] + FaceJ[i + 1][j].FK[4] - FaceJ[i][j].FK[4];
			//cout << Cell[i][j].QV[1] << " " << Cell[i][j].QV[2] << endl;
		}
	}
}
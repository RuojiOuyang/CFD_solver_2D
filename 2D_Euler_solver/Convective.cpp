/*
Convective函数计算各中心点的通量积分Q
*/

#include "Main.h"

using namespace std;

void Convective()
{
	//1.所有边界的通量
	//1.1边界处的通量
	//在Boundary_condition中给出
	
	//1.2内部边界处的通量

	//1.2.1内部边界处的守恒量
	for (j = 2; j <= Jmax - 1; j++)//i=1:Imax-1,j=2:Jmax-1对应的边都在内部
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			FaceI[i][j].FU[1] = 0.5 * (Cell[i][j - 1].U[1] + Cell[i][j].U[1]);//边界上的等于两边中心点求平均
			FaceI[i][j].FU[2] = 0.5 * (Cell[i][j - 1].U[2] + Cell[i][j].U[2]);
			FaceI[i][j].FU[3] = 0.5 * (Cell[i][j - 1].U[3] + Cell[i][j].U[3]);
			FaceI[i][j].FU[4] = 0.5 * (Cell[i][j - 1].U[4] + Cell[i][j].U[4]);
		}
	}
	for (j = 1; j <= Jmax - 1; j++)//i=2:Imax-1,j=1:Jmax-1对应的边都在内部
	{
		for (i = 2; i <= Imax - 1; i++)
		{
			FaceJ[i][j].FU[1] = 0.5 * (Cell[i - 1][j].U[1] + Cell[i][j].U[1]);//边界上的等于两边中心点求平均
			FaceJ[i][j].FU[2] = 0.5 * (Cell[i - 1][j].U[2] + Cell[i][j].U[2]);
			FaceJ[i][j].FU[3] = 0.5 * (Cell[i - 1][j].U[3] + Cell[i][j].U[3]);
			FaceJ[i][j].FU[4] = 0.5 * (Cell[i - 1][j].U[4] + Cell[i][j].U[4]);
		}
	}

	//1.2.2由内部边界上的守恒量计算内部边界上的通量
	double A, B;

	for (j = 2; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			A = FaceI[i][j].FU[1] * (pow(FaceI[i][j].FU[2] / FaceI[i][j].FU[1], 2.0) + pow(FaceI[i][j].FU[3] / FaceI[i][j].FU[1], 2.0)) / 2.0; //A是密度*动能
			B = (k - 1.0) * (FaceI[i][j].FU[4] - A);//B是压力
			FaceI[i][j].FF[1] = FaceI[i][j].FU[1] * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny);
			FaceI[i][j].FF[2] = FaceI[i][j].FU[2] * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny) + B * FaceI[i][j].nx;
			FaceI[i][j].FF[3] = FaceI[i][j].FU[3] * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny) + B * FaceI[i][j].ny;
			FaceI[i][j].FF[4] = (FaceI[i][j].FU[4] + B) * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny);
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 2; i <= Imax - 1; i++)
		{
			A = FaceJ[i][j].FU[1] * (pow(FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1], 2.0) + pow(FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1], 2.0)) / 2.0;
			B = (k - 1.0) * (FaceJ[i][j].FU[4] - A);
			FaceJ[i][j].FF[1] = FaceJ[i][j].FU[1] * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny);
			FaceJ[i][j].FF[2] = FaceJ[i][j].FU[2] * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny) + B * FaceJ[i][j].nx;
			FaceJ[i][j].FF[3] = FaceJ[i][j].FU[3] * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny) + B * FaceJ[i][j].ny;
			FaceJ[i][j].FF[4] = (FaceJ[i][j].FU[4] + B) * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny);
		}
	}

	//2.有了所有边界的通量后，控制体单元的通量，也就是Qij
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].QC[1] = FaceI[i][j + 1].FF[1] - FaceI[i][j].FF[1] + FaceJ[i + 1][j].FF[1] - FaceJ[i][j].FF[1];
			Cell[i][j].QC[2] = FaceI[i][j + 1].FF[2] - FaceI[i][j].FF[2] + FaceJ[i + 1][j].FF[2] - FaceJ[i][j].FF[2];
			Cell[i][j].QC[3] = FaceI[i][j + 1].FF[3] - FaceI[i][j].FF[3] + FaceJ[i + 1][j].FF[3] - FaceJ[i][j].FF[3];
			Cell[i][j].QC[4] = FaceI[i][j + 1].FF[4] - FaceI[i][j].FF[4] + FaceJ[i + 1][j].FF[4] - FaceJ[i][j].FF[4];
		}
	}
}
/*
Local_time_step函数计算当地时间步
*/

#include "Main.h"

using namespace std;

void Min_time_step()
{
	double A, B, C, D, E, F, G, L;//计算当地时间步时的中间变量，局部变量
	double mintime;

	mintime = 100;

	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			A = 0.5 * (FaceJ[i][j].nx + FaceJ[i + 1][j].nx);//这个式子就是把式(2-1-13)里面第一项的SI里的x方向的加在一起了
			B = 0.5 * (FaceJ[i][j].ny + FaceJ[i + 1][j].ny);//这个式子就是把式(2-1-13)里面第一项的SI里的y方向的加在一起了
			C = 0.5 * (FaceI[i][j].nx + FaceI[i][j + 1].nx);//这个式子就是把式(2-1-13)里面第二项的SI里的x方向的加在一起了
			D = 0.5 * (FaceI[i][j].ny + FaceI[i][j + 1].ny);//这个式子就是把式(2-1-13)里面第二项的SI里的y方向的加在一起了
			E = fabs(Cell[i][j].u * A + Cell[i][j].v * B);//式(2-1-13)里面第一项
			F = fabs(Cell[i][j].u * C + Cell[i][j].v * D);//式(2-1-13)里面第二项
			G = fabs(sqrt(A * A + B * B));
			L = fabs(sqrt(C * C + D * D));//式(2-1-13)里面第三项
			Cell[i][j].localdt = CFL * Cell[i][j].Vol / (E + F + Cell[i][j].c * (G + L));//计算出当地时间步，存在中心点处

			if (Cell[i][j].localdt < mintime)//记录下当前最小时间步
			{
				mintime = Cell[i][j].localdt;
			}
		}
	}
	
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].dt = mintime;//推进时间采用最小时间步
		}
	}

	totaltime += mintime;
}
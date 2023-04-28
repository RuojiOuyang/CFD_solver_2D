/*
Solve_conservation函数利用变量求通量，求流场上各点（中心点）的通量
*/

#include "Main.h"

using namespace std;

void Solve_conservation()
{
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].U[1] = Cell[i][j].Den;
			Cell[i][j].U[2] = Cell[i][j].Den * Cell[i][j].u;
			Cell[i][j].U[3] = Cell[i][j].Den * Cell[i][j].v;
			Cell[i][j].U[4] = Cell[i][j].Den * Cell[i][j].E;
		}
	}
}
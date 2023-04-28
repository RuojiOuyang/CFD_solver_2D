/*
Primitive函数从守恒量反求原始变量，用于下步计算赋值
*/

#include "Main.h"

using namespace std;

void Primitive()
{
	double Kin;//动能

	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].Den = Cell[i][j].U[1];
			Cell[i][j].u = Cell[i][j].U[2] / Cell[i][j].U[1];
			Cell[i][j].v = Cell[i][j].U[3] / Cell[i][j].U[1];
			Cell[i][j].E = Cell[i][j].U[4] / Cell[i][j].U[1];
			Kin = (Cell[i][j].u * Cell[i][j].u + Cell[i][j].v * Cell[i][j].v) / 2.0;
			Cell[i][j].P = (k - 1.0) * Cell[i][j].Den * (Cell[i][j].E - Kin);
			Cell[i][j].H = Cell[i][j].E + Cell[i][j].P / Cell[i][j].Den;
			Cell[i][j].T = Cell[i][j].P * k / Cell[i][j].Den;
			//if ((Cell[i][j].T) < 0)
			//{
			//	cout << "error" << endl;
			//}
			Cell[i][j].c = sqrt(Cell[i][j].T);
			Cell[i][j].Ma = sqrt(Kin * 2.0) / Cell[i][j].c;
		}
	}
}
/*
Solve_conservation�������ñ�����ͨ�����������ϸ��㣨���ĵ㣩��ͨ��
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
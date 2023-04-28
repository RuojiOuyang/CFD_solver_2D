/*
Read_mesh函数用来读网格，把网格点数据存在Node里面
*/

#include "Main.h"

using namespace std;

void ReadMesh()
{
	ifstream read;
	read.open("naca0012_09_17data.txt");

	read >> Imax >> Jmax;//读前两个数，前两个数分别是i,j方向网格点的数量

	for (j = 0; j < Jmax; j++)
	{
		for (i = 0; i < Imax; i++)
		{
			read >> Node_x(i, j) >> Node_y(i, j);
		}
	}
	read.close();
}
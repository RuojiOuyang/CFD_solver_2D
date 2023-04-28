/*
Read_mesh函数用来读网格，把网格点数据存在Node里面
*/

#include "Main.h"

using namespace std;

void Read_mesh()
{
	ifstream read;
	read.open("yuanzhudata.txt");

	read >> Imax >> Jmax;//读前两个数，前两个数分别是i,j方向网格点的数量

	for (j = 1; j <= Jmax; j++)
	{
		for (i = 1; i <= Imax; i++)
		{
			read >> Node[i][j].x >> Node[i][j].y;
			Node[i][j].x = Node[i][j].x * 1e-3;
			Node[i][j].y = Node[i][j].y * 1e-3;
			/*if (Node[i][j].x > 1 || Node[i][j].y > 1)
			{
				exit(0);
			}*/
		}
	}
	read.close();
}
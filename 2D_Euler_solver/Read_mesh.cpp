/*
Read_mesh�������������񣬰���������ݴ���Node����
*/

#include "Main.h"

using namespace std;

void Read_mesh()
{
	ifstream read;
	read.open("naca0012data.txt");

	read >> Imax >> Jmax;//��ǰ��������ǰ�������ֱ���i,j��������������

	for (j = 1; j <= Jmax; j++)
	{
		for (i = 1; i <= Imax; i++)
		{
			read >> Node[i][j].x >> Node[i][j].y;
			Node[i][j].x = Node[i][j].x * 1e-1;
			Node[i][j].y = Node[i][j].y * 1e-1;
			if (Node[i][j].x > 100 || Node[i][j].y > 100)
			{
				exit(0);
			}
		}
	}
	read.close();
}
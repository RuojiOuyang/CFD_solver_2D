/*
Read_mesh�������������񣬰���������ݴ���Node����
*/

#include "Main.h"

using namespace std;

void ReadMesh()
{
	ifstream read;
	read.open("naca0012_09_17data.txt");

	read >> Imax >> Jmax;//��ǰ��������ǰ�������ֱ���i,j��������������

	for (j = 0; j < Jmax; j++)
	{
		for (i = 0; i < Imax; i++)
		{
			read >> Node_x(i, j) >> Node_y(i, j);
		}
	}
	read.close();
}
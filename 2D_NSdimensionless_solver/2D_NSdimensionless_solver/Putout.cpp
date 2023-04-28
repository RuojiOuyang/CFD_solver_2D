/*
Putout函数，将结果输出至文件
*/

#include "Main.h"

using namespace std;

void Putout()
{
	ofstream out;
	out.open("result.dat");
	out << "variables = x y Den P T Ma u v E H" << endl;
	out << "zone i=" << Imax << " and j=" << Jmax << endl;
	for (j = 1; j <= Jmax; j++)
	{
		for (i = 1; i <= Imax; i++)
		{
			out << Node[i][j].x * L << " " << Node[i][j].y * L << " " << Node[i][j].Den * Dentinf << " " << Node[i][j].P * Dentinf * ctinf * ctinf << " " << Node[i][j].T * Ttinf << " " << Node[i][j].Ma << " " << Node[i][j].u * ctinf << " " << Node[i][j].v * vtinf << " " << Node[i][j].E * ctinf * ctinf << " " << Node[i][j].H * ctinf * ctinf << " " << Node[i][j].u_y * ctinf / L << " " << Node[i][j].v_x * ctinf / L << endl;
		}
	}
	out.close();
}
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
			out.precision(16);
			out << Node[i][j].x << " " << Node[i][j].y << " " << Node[i][j].Den << " " << Node[i][j].P << " " << Node[i][j].T << " " << Node[i][j].Ma << " " << Node[i][j].u << " " << Node[i][j].v << " " << Node[i][j].E << " " << Node[i][j].H << " " << Node[i][j].u_y << " " << Node[i][j].v_x << endl;
		}
	}
	out.close();
}
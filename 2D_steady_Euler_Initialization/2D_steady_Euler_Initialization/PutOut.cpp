/*
PutOut函数，将结果输出至文件
*/

#include "Main.h"

using namespace std;
using namespace arma;

void PutOut()
{
	string s = "SteadyInitialization";
	string s1 = ".dat";
	stringstream ss;

	ss << s << numAOA << s1;

	ofstream out;
	out.open(ss.str());
	//out << "variables = x y Den P T Ma u v E" << endl;
	//out << "zone i=" << Imax << " and j=" << Jmax << endl;
	/*for (j = 0; j < Jmax; j++)
	{
		for (i = 0; i < Imax; i++)
		{
			out << Node_x(i, j) << " " << Node_y(i, j) << " "
				<< Node_Den(i, j) << " " << Node_P(i, j) << " "
				<< Node_T(i, j) << " " << Node_Ma(i, j) << " "
				<< Node_u(i, j) << " " << Node_v(i, j) << " "
				<< Node_E(i, j) << endl;
		}
	}*/
	for (j = 0; j < Jmax - 1; j++)
	{
		for (i = 0; i < Imax - 1; i++)
		{
			out << Cell_Den(i, j) << " " << Cell_P(i, j) << " "
				<< Cell_T(i, j) << " " << Cell_Ma(i, j) << " "
				<< Cell_u(i, j) << " " << Cell_v(i, j) << " "
				<< Cell_E(i, j) << endl;
		}
	}
	out.close();

	ss.clear();//恢复状态
	ss.str("");//恢复值
}
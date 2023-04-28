/*
Primitive函数从守恒量反求原始变量，用于下步计算赋值
*/

#include "Main.h"

using namespace std;
using namespace arma;

void Primitive()
{
	mat Kin;//动能

	Cell_Den = Cell_U1;
	Cell_u = Cell_U2 / Cell_U1;
	Cell_v = Cell_U3 / Cell_U1;
	Cell_E = Cell_U4 / Cell_U1;
	Kin = (Cell_u % Cell_u + Cell_v % Cell_v) / 2.0;
	Cell_P = (k - 1.0) * Cell_Den % (Cell_E - Kin);
	Cell_T = Cell_P / R / Cell_Den;
	Cell_c = arma::sqrt(k * R * Cell_T);
	Cell_Ma = arma::sqrt(Kin * 2.0) / Cell_c;
}
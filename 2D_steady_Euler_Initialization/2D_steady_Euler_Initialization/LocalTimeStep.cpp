/*
LocalTimeStep函数计算当地时间步，参数为时间层
*/

#include "Main.h"

using namespace std;
using namespace arma;

void LocalTimeStep()
{
	//计算当地时间步时的中间变量
	mat AA, BB, CC, DD;
	mat EE, FF, GG, LL, OO, PP, RR;

	AA = 0.5 * (FaceJ_nx.rows(0, (long long)Imax - 2) + FaceJ_nx.rows(1, (long long)Imax - 1));//这个式子就是把式(2-1-13)里面第一项的SI里的x方向的加在一起了
	BB = 0.5 * (FaceJ_ny.rows(0, (long long)Imax - 2) + FaceJ_ny.rows(1, (long long)Imax - 1));//这个式子就是把式(2-1-13)里面第一项的SI里的y方向的加在一起了
	CC = 0.5 * (FaceI_nx.cols(0, (long long)Jmax - 2) + FaceI_nx.cols(1, (long long)Jmax - 1));//这个式子就是把式(2-1-13)里面第二项的SI里的x方向的加在一起了
	DD = 0.5 * (FaceI_ny.cols(0, (long long)Jmax - 2) + FaceI_ny.cols(1, (long long)Jmax - 1));//这个式子就是把式(2-1-13)里面第二项的SI里的y方向的加在一起了

	EE = arma::abs(Cell_u % AA + Cell_v % BB);//式(2-1-13)里面第一项
	FF = arma::abs(Cell_u % CC + Cell_v % DD);//式(2-1-13)里面第二项
	GG = arma::abs(arma::sqrt(AA % AA + BB % BB));
	LL = arma::abs(arma::sqrt(CC % CC + DD % DD));
	OO = EE + Cell_c % GG;
	PP = FF + Cell_c % LL;

	RR = OO + PP;

	Cell_dt = CFL * Cell_Vol / RR;//计算出当地时间步，存在中心点处
}
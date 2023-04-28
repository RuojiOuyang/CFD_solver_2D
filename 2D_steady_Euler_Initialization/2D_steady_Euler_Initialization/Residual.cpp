/*
Residualº¯Êý¼ÆËã²Ð²î
*/

#include "Main.h"

using namespace std;
using namespace arma;

void Residual()
{
	Resid = 0.0;

	Resid = Resid + accu(arma::pow(P_former - Cell_P, 2.0));

	Resid = sqrt(Resid / (Imax - 1.0) / (Jmax - 1.0));
}
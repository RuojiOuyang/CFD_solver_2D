/*
Residualº¯Êý¼ÆËã²Ð²î
*/

#include "Main.h"

using namespace std;

void Residual()
{
	Resid = 0.0;
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			Resid = Resid + pow(Den_former[i][j] - Cell[i][j].Den, 2.0);
		}
	}
	Resid = sqrt(Resid / (Imax - 1.0) / (Jmax - 1.0));
}
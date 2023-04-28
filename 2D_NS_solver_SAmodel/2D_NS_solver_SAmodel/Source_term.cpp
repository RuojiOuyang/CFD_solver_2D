/*
Source_term项计算SA湍流模型产生的源项
*/

#include "Main.h"

using namespace std;

void Source_term()
{
	//要计算源项的第三项，首先需要计算边界上的偏导数的值，网格边界上的用虚元与内部中心差分得到
	//1.偏导数值
	//1.1边界边界
	//1.1.1下上边界
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][1].miubl_x = 0.5 * (Cell[i][1].miubl_x + DummyIS[i][1].miubl_x);
		FaceI[i][1].miubl_y = 0.5 * (Cell[i][1].miubl_y + DummyIS[i][1].miubl_y);
		FaceI[i][Jmax].miubl_x = 0.5 * (Cell[i][Jmax - 1].miubl_x + DummyIN[i][1].miubl_x);
		FaceI[i][Jmax].miubl_y = 0.5 * (Cell[i][Jmax - 1].miubl_y + DummyIN[i][1].miubl_y);
	}
	//1.1.2左右边界
	for (j = 1; j <= Jmax - 1; j++)
	{
		FaceJ[1][j].miubl_x = 0.5 * (Cell[1][j].miubl_x + DummyJW[1][j].miubl_x);
		FaceJ[1][j].miubl_y = 0.5 * (Cell[1][j].miubl_y + DummyJW[1][j].miubl_y);
		FaceJ[Imax][j].miubl_x = 0.5 * (Cell[Imax - 1][j].miubl_x + DummyJE[1][j].miubl_x);
		FaceJ[Imax][j].miubl_y = 0.5 * (Cell[Imax - 1][j].miubl_y + DummyJE[1][j].miubl_y);
	}
	//1.2内部边界
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 2; j <= Jmax - 1; j++)
		{
			FaceI[i][j].miubl_x = 0.5 * (Cell[i][j].miubl_x + Cell[i][j - 1].miubl_x);
			FaceI[i][j].miubl_y = 0.5 * (Cell[i][j].miubl_y + Cell[i][j - 1].miubl_y);
		}
	}
	for (i = 2; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			FaceJ[i][j].miubl_x = 0.5 * (Cell[i][j].miubl_x + Cell[i - 1][j].miubl_x);
			FaceJ[i][j].miubl_y = 0.5 * (Cell[i][j].miubl_y + Cell[i - 1][j].miubl_y);
		}
	}

	//2源项的计算
	double RESULT;
	double ft2;
	double fv2;
	double S;
	double Sbl;//S一波浪
	double Ct3 = 1.2;
	double Ct4 = 0.5;
	double fv3 = 1.0;
	//double Cv2 = 5;
	double fw, g, r;

	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			RESULT = 0;
			
			ft2 = Ct3 * exp(-Ct4 * Cell[i][j].writeX * Cell[i][j].writeX);
			fv2 = 1 - Cell[i][j].writeX / (1 + Cell[i][j].writeX * Cell[i][j].fv1);
			S = sqrt(2 * (1.0 / 2.0 * (Cell[i][j].u_y - Cell[i][j].v_x)) * (1.0 / 2.0 * (Cell[i][j].u_y - Cell[i][j].v_x))) * fv3;
			Sbl = S + Cell[i][j].U[5] / Cell[i][j].U[1] / writek / writek / Cell[i][j].d / Cell[i][j].d * fv2;
			Cell[i][j].S = Cb1 * ((1 - ft2) * Sbl * Cell[i][j].U[5]);//第一项

			r = min(Cell[i][j].U[5] / Cell[i][j].U[1] / writek / writek / Cell[i][j].d / Cell[i][j].d / Sbl, rlim);
			g = r + Cw2 * (pow(r, 6) - r);
			fw = g * pow((1 + pow(Cw3, 6)) / (pow(g, 6) + pow(Cw3, 6)), 1.0 / 6.0);
			Cell[i][j].S -= (Cw1 * fw - Cb1 / writek / writek * ft2) * Cell[i][j].U[1] * Cell[i][j].U[5] / Cell[i][j].U[1] * Cell[i][j].U[5] / Cell[i][j].U[1] / Cell[i][j].d / Cell[i][j].d;//第二项
			Cell[i][j].S = Cell[i][j].S * Cell[i][j].Vol;

			RESULT = FaceI[i][j + 1].miubl_x * FaceI[i][j + 1].nx + FaceI[i][j + 1].miubl_y * FaceI[i][j + 1].ny;
			RESULT += FaceJ[i + 1][j].miubl_x * FaceJ[i + 1][j].nx + FaceJ[i + 1][j].miubl_y * FaceJ[i + 1][j].ny;
			RESULT -= FaceI[i][j].miubl_x * FaceI[i][j].nx + FaceI[i][j].miubl_y * FaceI[i][j].ny;
			RESULT -= FaceJ[i][j].miubl_x * FaceJ[i][j].nx + FaceJ[i][j].miubl_y * FaceJ[i][j].ny;
			Cell[i][j].S -= Cb2 / delta * Cell[i][j].U[5] * RESULT;//第三项
		}
	}
}
/*
Calculate_P函数计算T、u、v在xy方向的偏导数
*/

#include "Main.h"

using namespace std;

void Calculate_P()
{
	double A1, A2, A3, A4;
	double B1, B2, B3, B4;//计算过程中的中间量，局部变量

	//1.四个角点
	{
		//1.1T
		{
			FaceI[1][1].T = 0.5 * (Cell[1][1].T + DummyIS[1][1].T);
			FaceJ[1][1].T = 0.5 * (Cell[1][1].T + DummyJW[1][1].T);
			FaceI[Imax - 1][1].T = 0.5 * (Cell[Imax - 1][1].T + DummyIS[Imax - 1][1].T);
			FaceJ[Imax][1].T = 0.5 * (Cell[Imax - 1][1].T + DummyJE[1][1].T);
			FaceI[Imax - 1][Jmax].T = 0.5 * (Cell[Imax - 1][Jmax - 1].T + DummyIN[Imax - 1][1].T);
			FaceJ[Imax][Jmax - 1].T = 0.5 * (Cell[Imax - 1][Jmax - 1].T + DummyJE[1][Jmax - 1].T);
			FaceI[1][Jmax].T = 0.5 * (Cell[1][Jmax - 1].T + DummyIN[1][1].T);
			FaceJ[1][Jmax - 1].T = 0.5 * (Cell[1][Jmax - 1].T + DummyJW[1][Jmax - 1].T);

			A1 = (Cell[1][1].T + Cell[2][1].T) * FaceJ[2][1].nx / 2;
			A2 = FaceJ[1][1].T * FaceJ[1][1].nx;
			A3 = (Cell[1][1].T + Cell[1][2].T) * FaceI[1][2].nx / 2;
			A4 = FaceI[1][1].T * FaceI[1][1].nx;
			Cell[1][1].T_x = (A1 - A2 + A3 - A4) / Cell[1][1].Vol;
			B1 = (Cell[1][1].T + Cell[2][1].T) * FaceJ[2][1].ny / 2;
			B2 = FaceJ[1][1].T * FaceJ[1][1].ny;
			B3 = (Cell[1][1].T + Cell[1][2].T) * FaceI[1][2].ny / 2;
			B4 = FaceI[1][1].T * FaceI[1][1].ny;
			Cell[1][1].T_y = (B1 - B2 + B3 - B4) / Cell[1][1].Vol;
			A1 = (Cell[1][Jmax - 1].T + Cell[2][Jmax - 1].T) * FaceJ[2][Jmax - 1].nx / 2;
			A2 = FaceJ[1][Jmax - 1].T * FaceJ[1][Jmax - 1].nx;
			A3 = FaceI[1][Jmax].T * FaceI[1][Jmax].nx;
			A4 = (Cell[1][Jmax - 1].T + Cell[1][Jmax - 2].T) * FaceI[1][Jmax - 1].nx / 2;
			Cell[1][Jmax - 1].T_x = (A1 - A2 + A3 - A4) / Cell[1][Jmax - 1].Vol;
			B1 = (Cell[1][Jmax - 1].T + Cell[2][Jmax - 1].T) * FaceJ[2][Jmax - 1].ny / 2;
			B2 = FaceJ[1][Jmax - 1].T * FaceJ[1][Jmax - 1].ny;
			B3 = FaceI[1][Jmax].T * FaceI[1][Jmax].ny;
			B4 = (Cell[1][Jmax - 1].T + Cell[1][Jmax - 2].T) * FaceI[1][Jmax - 1].ny / 2;
			Cell[1][Jmax - 1].T_y = (B1 - B2 + B3 - B4) / Cell[1][Jmax - 1].Vol;
			A1 = FaceJ[Imax][1].T * FaceJ[Imax][1].nx;
			A2 = (Cell[Imax - 1][1].T + Cell[Imax - 2][1].T) * FaceJ[Imax - 1][1].nx / 2;
			A3 = (Cell[Imax - 1][1].T + Cell[Imax - 1][2].T) * FaceI[Imax - 1][2].nx / 2;
			A4 = FaceI[Imax - 1][1].T * FaceI[Imax - 1][1].nx;
			Cell[Imax - 1][1].T_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][1].Vol;
			B1 = FaceJ[Imax][1].T * FaceJ[Imax][1].ny;
			B2 = (Cell[Imax - 1][1].T + Cell[Imax - 2][1].T) * FaceJ[Imax - 1][1].ny / 2;
			B3 = (Cell[Imax - 1][1].T + Cell[Imax - 1][2].T) * FaceI[Imax - 1][2].ny / 2;
			B4 = FaceI[Imax - 1][1].T * FaceI[Imax - 1][1].ny;
			Cell[Imax - 1][1].T_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][1].Vol;
			A1 = FaceJ[Imax][Jmax - 1].T * FaceJ[Imax][Jmax - 1].nx;
			A2 = (Cell[Imax - 1][Jmax - 1].T + Cell[Imax - 2][Jmax - 1].T) * FaceJ[Imax - 1][Jmax - 1].nx / 2;
			A3 = FaceI[Imax - 1][Jmax].T * FaceI[Imax - 1][Jmax].nx;
			A4 = (Cell[Imax - 1][Jmax - 1].T + Cell[Imax - 1][Jmax - 2].T) * FaceI[Imax - 1][Jmax - 1].nx / 2;
			Cell[Imax - 1][Jmax - 1].T_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][Jmax - 1].Vol;
			B1 = FaceJ[Imax][Jmax - 1].T * FaceJ[Imax][Jmax - 1].ny;
			B2 = (Cell[Imax - 1][Jmax - 1].T + Cell[Imax - 2][Jmax - 1].T) * FaceJ[Imax - 1][Jmax - 1].ny / 2;
			B3 = FaceI[Imax - 1][Jmax].T * FaceI[Imax - 1][Jmax].ny;
			B4 = (Cell[Imax - 1][Jmax - 1].T + Cell[Imax - 1][Jmax - 2].T) * FaceI[Imax - 1][Jmax - 1].ny / 2;
			Cell[Imax - 1][Jmax - 1].T_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][Jmax - 1].Vol;
		}
		//1.2u
		{
			FaceI[1][1].u = 0.5 * (Cell[1][1].u + DummyIS[1][1].u);
			FaceJ[1][1].u = 0.5 * (Cell[1][1].u + DummyJW[1][1].u);
			FaceI[Imax - 1][1].u = 0.5 * (Cell[Imax - 1][1].u + DummyIS[Imax - 1][1].u);
			FaceJ[Imax][1].u = 0.5 * (Cell[Imax - 1][1].u + DummyJE[1][1].u);
			FaceI[Imax - 1][Jmax].u = 0.5 * (Cell[Imax - 1][Jmax - 1].u + DummyIN[Imax - 1][1].u);
			FaceJ[Imax][Jmax - 1].u = 0.5 * (Cell[Imax - 1][Jmax - 1].u + DummyJE[1][Jmax - 1].u);
			FaceI[1][Jmax].u = 0.5 * (Cell[1][Jmax - 1].u + DummyIN[1][1].u);
			FaceJ[1][Jmax - 1].u = 0.5 * (Cell[1][Jmax - 1].u + DummyJW[1][Jmax - 1].u);

			A1 = (Cell[1][1].u + Cell[2][1].u) * FaceJ[2][1].nx / 2;
			A2 = FaceJ[1][1].u * FaceJ[1][1].nx;
			A3 = (Cell[1][1].u + Cell[1][2].u) * FaceI[1][2].nx / 2;
			A4 = FaceI[1][1].u * FaceI[1][1].nx;
			Cell[1][1].u_x = (A1 - A2 + A3 - A4) / Cell[1][1].Vol;
			B1 = (Cell[1][1].u + Cell[2][1].u) * FaceJ[2][1].ny / 2;
			B2 = FaceJ[1][1].u * FaceJ[1][1].ny;
			B3 = (Cell[1][1].u + Cell[1][2].u) * FaceI[1][2].ny / 2;
			B4 = FaceI[1][1].u * FaceI[1][1].ny;
			Cell[1][1].u_y = (B1 - B2 + B3 - B4) / Cell[1][1].Vol;
			A1 = (Cell[1][Jmax - 1].u + Cell[2][Jmax - 1].u) * FaceJ[2][Jmax - 1].nx / 2;
			A2 = FaceJ[1][Jmax - 1].u * FaceJ[1][Jmax - 1].nx;
			A3 = FaceI[1][Jmax].u * FaceI[1][Jmax].nx;
			A4 = (Cell[1][Jmax - 1].u + Cell[1][Jmax - 2].u) * FaceI[1][Jmax - 1].nx / 2;
			Cell[1][Jmax - 1].u_x = (A1 - A2 + A3 - A4) / Cell[1][Jmax - 1].Vol;
			B1 = (Cell[1][Jmax - 1].u + Cell[2][Jmax - 1].u) * FaceJ[2][Jmax - 1].ny / 2;
			B2 = FaceJ[1][Jmax - 1].u * FaceJ[1][Jmax - 1].ny;
			B3 = FaceI[1][Jmax].u * FaceI[1][Jmax].ny;
			B4 = (Cell[1][Jmax - 1].u + Cell[1][Jmax - 2].u) * FaceI[1][Jmax - 1].ny / 2;
			Cell[1][Jmax - 1].u_y = (B1 - B2 + B3 - B4) / Cell[1][Jmax - 1].Vol;
			A1 = FaceJ[Imax][1].u * FaceJ[Imax][1].nx;
			A2 = (Cell[Imax - 1][1].u + Cell[Imax - 2][1].u) * FaceJ[Imax - 1][1].nx / 2;
			A3 = (Cell[Imax - 1][1].u + Cell[Imax - 1][2].u) * FaceI[Imax - 1][2].nx / 2;
			A4 = FaceI[Imax - 1][1].u * FaceI[Imax - 1][1].nx;
			Cell[Imax - 1][1].u_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][1].Vol;
			B1 = FaceJ[Imax][1].u * FaceJ[Imax][1].ny;
			B2 = (Cell[Imax - 1][1].u + Cell[Imax - 2][1].u) * FaceJ[Imax - 1][1].ny / 2;
			B3 = (Cell[Imax - 1][1].u + Cell[Imax - 1][2].u) * FaceI[Imax - 1][2].ny / 2;
			B4 = FaceI[Imax - 1][1].u * FaceI[Imax - 1][1].ny;
			Cell[Imax - 1][1].u_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][1].Vol;
			A1 = FaceJ[Imax][Jmax - 1].u * FaceJ[Imax][Jmax - 1].nx;
			A2 = (Cell[Imax - 1][Jmax - 1].u + Cell[Imax - 2][Jmax - 1].u) * FaceJ[Imax - 1][Jmax - 1].nx / 2;
			A3 = FaceI[Imax - 1][Jmax].u * FaceI[Imax - 1][Jmax].nx;
			A4 = (Cell[Imax - 1][Jmax - 1].u + Cell[Imax - 1][Jmax - 2].u) * FaceI[Imax - 1][Jmax - 1].nx / 2;
			Cell[Imax - 1][Jmax - 1].u_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][Jmax - 1].Vol;
			B1 = FaceJ[Imax][Jmax - 1].u * FaceJ[Imax][Jmax - 1].ny;
			B2 = (Cell[Imax - 1][Jmax - 1].u + Cell[Imax - 2][Jmax - 1].u) * FaceJ[Imax - 1][Jmax - 1].ny / 2;
			B3 = FaceI[Imax - 1][Jmax].u * FaceI[Imax - 1][Jmax].ny;
			B4 = (Cell[Imax - 1][Jmax - 1].u + Cell[Imax - 1][Jmax - 2].u) * FaceI[Imax - 1][Jmax - 1].ny / 2;
			Cell[Imax - 1][Jmax - 1].u_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][Jmax - 1].Vol;
		}
		//1.3v
		{
			FaceI[1][1].v = 0.5 * (Cell[1][1].v + DummyIS[1][1].v);
			FaceJ[1][1].v = 0.5 * (Cell[1][1].v + DummyJW[1][1].v);
			FaceI[Imax - 1][1].v = 0.5 * (Cell[Imax - 1][1].v + DummyIS[Imax - 1][1].v);
			FaceJ[Imax][1].v = 0.5 * (Cell[Imax - 1][1].v + DummyJE[1][1].v);
			FaceI[Imax - 1][Jmax].v = 0.5 * (Cell[Imax - 1][Jmax - 1].v + DummyIN[Imax - 1][1].v);
			FaceJ[Imax][Jmax - 1].v = 0.5 * (Cell[Imax - 1][Jmax - 1].v + DummyJE[1][Jmax - 1].v);
			FaceI[1][Jmax].v = 0.5 * (Cell[1][Jmax - 1].v + DummyIN[1][1].v);
			FaceJ[1][Jmax - 1].v = 0.5 * (Cell[1][Jmax - 1].v + DummyJW[1][Jmax - 1].v);

			A1 = (Cell[1][1].v + Cell[2][1].v) * FaceJ[2][1].nx / 2;
			A2 = FaceJ[1][1].v * FaceJ[1][1].nx;
			A3 = (Cell[1][1].v + Cell[1][2].v) * FaceI[1][2].nx / 2;
			A4 = FaceI[1][1].v * FaceI[1][1].nx;
			Cell[1][1].v_x = (A1 - A2 + A3 - A4) / Cell[1][1].Vol;
			B1 = (Cell[1][1].v + Cell[2][1].v) * FaceJ[2][1].ny / 2;
			B2 = FaceJ[1][1].v * FaceJ[1][1].ny;
			B3 = (Cell[1][1].v + Cell[1][2].v) * FaceI[1][2].ny / 2;
			B4 = FaceI[1][1].v * FaceI[1][1].ny;
			Cell[1][1].v_y = (B1 - B2 + B3 - B4) / Cell[1][1].Vol;
			A1 = (Cell[1][Jmax - 1].v + Cell[2][Jmax - 1].v) * FaceJ[2][Jmax - 1].nx / 2;
			A2 = FaceJ[1][Jmax - 1].v * FaceJ[1][Jmax - 1].nx;
			A3 = FaceI[1][Jmax].v * FaceI[1][Jmax].nx;
			A4 = (Cell[1][Jmax - 1].v + Cell[1][Jmax - 2].v) * FaceI[1][Jmax - 1].nx / 2;
			Cell[1][Jmax - 1].v_x = (A1 - A2 + A3 - A4) / Cell[1][Jmax - 1].Vol;
			B1 = (Cell[1][Jmax - 1].v + Cell[2][Jmax - 1].v) * FaceJ[2][Jmax - 1].ny / 2;
			B2 = FaceJ[1][Jmax - 1].v * FaceJ[1][Jmax - 1].ny;
			B3 = FaceI[1][Jmax].v * FaceI[1][Jmax].ny;
			B4 = (Cell[1][Jmax - 1].v + Cell[1][Jmax - 2].v) * FaceI[1][Jmax - 1].ny / 2;
			Cell[1][Jmax - 1].v_y = (B1 - B2 + B3 - B4) / Cell[1][Jmax - 1].Vol;
			A1 = FaceJ[Imax][1].v * FaceJ[Imax][1].nx;
			A2 = (Cell[Imax - 1][1].v + Cell[Imax - 2][1].v) * FaceJ[Imax - 1][1].nx / 2;
			A3 = (Cell[Imax - 1][1].v + Cell[Imax - 1][2].v) * FaceI[Imax - 1][2].nx / 2;
			A4 = FaceI[Imax - 1][1].v * FaceI[Imax - 1][1].nx;
			Cell[Imax - 1][1].v_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][1].Vol;
			B1 = FaceJ[Imax][1].v * FaceJ[Imax][1].ny;
			B2 = (Cell[Imax - 1][1].v + Cell[Imax - 2][1].v) * FaceJ[Imax - 1][1].ny / 2;
			B3 = (Cell[Imax - 1][1].v + Cell[Imax - 1][2].v) * FaceI[Imax - 1][2].ny / 2;
			B4 = FaceI[Imax - 1][1].v * FaceI[Imax - 1][1].ny;
			Cell[Imax - 1][1].v_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][1].Vol;
			A1 = FaceJ[Imax][Jmax - 1].v * FaceJ[Imax][Jmax - 1].nx;
			A2 = (Cell[Imax - 1][Jmax - 1].v + Cell[Imax - 2][Jmax - 1].v) * FaceJ[Imax - 1][Jmax - 1].nx / 2;
			A3 = FaceI[Imax - 1][Jmax].v * FaceI[Imax - 1][Jmax].nx;
			A4 = (Cell[Imax - 1][Jmax - 1].v + Cell[Imax - 1][Jmax - 2].v) * FaceI[Imax - 1][Jmax - 1].nx / 2;
			Cell[Imax - 1][Jmax - 1].v_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][Jmax - 1].Vol;
			B1 = FaceJ[Imax][Jmax - 1].v * FaceJ[Imax][Jmax - 1].ny;
			B2 = (Cell[Imax - 1][Jmax - 1].v + Cell[Imax - 2][Jmax - 1].v) * FaceJ[Imax - 1][Jmax - 1].ny / 2;
			B3 = FaceI[Imax - 1][Jmax].v * FaceI[Imax - 1][Jmax].ny;
			B4 = (Cell[Imax - 1][Jmax - 1].v + Cell[Imax - 1][Jmax - 2].v) * FaceI[Imax - 1][Jmax - 1].ny / 2;
			Cell[Imax - 1][Jmax - 1].v_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][Jmax - 1].Vol;
		}
	}

	//2.除四个角点外的四个边
	for (i = 2; i <= Imax - 2; i++)
	{
		FaceI[i][1].T = 0.5 * (Cell[i][1].T + DummyIS[i][1].T);
		FaceI[i][Jmax].T = 0.5 * (Cell[i][Jmax - 1].T + DummyIN[i][1].T);
		FaceI[i][1].u = 0.5 * (Cell[i][1].u + DummyIS[i][1].u);
		FaceI[i][Jmax].u = 0.5 * (Cell[i][Jmax - 1].u + DummyIN[i][1].u);
		FaceI[i][1].v = 0.5 * (Cell[i][1].v + DummyIS[i][1].v);
		FaceI[i][Jmax].v = 0.5 * (Cell[i][Jmax - 1].v + DummyIN[i][1].v);

		A1 = (Cell[i][1].T + Cell[i + 1][1].T) * FaceJ[i + 1][1].nx / 2;
		A2 = (Cell[i][1].T + Cell[i - 1][1].T) * FaceJ[i][1].nx / 2;
		A3 = (Cell[i][1].T + Cell[i][2].T) * FaceI[i][2].nx / 2;
		A4 = FaceI[i][1].T * FaceI[i][1].nx;
		Cell[i][1].T_x = (A1 - A2 + A3 - A4) / Cell[i][1].Vol;
		B1 = (Cell[i][1].T + Cell[i + 1][1].T) * FaceJ[i + 1][1].ny / 2;
		B2 = (Cell[i][1].T + Cell[i - 1][1].T) * FaceJ[i][1].ny / 2;
		B3 = (Cell[i][1].T + Cell[i][2].T) * FaceI[i][2].ny / 2;
		B4 = FaceI[i][1].T * FaceI[i][1].ny;
		Cell[i][1].T_y = (B1 - B2 + B3 - B4) / Cell[i][1].Vol;
		A1 = (Cell[i][Jmax - 1].T + Cell[i + 1][Jmax - 1].T) * FaceJ[i + 1][Jmax - 1].nx / 2;
		A2 = (Cell[i][Jmax - 1].T + Cell[i - 1][Jmax - 1].T) * FaceJ[i][Jmax - 1].nx / 2;
		A3 = FaceI[i][Jmax].T * FaceI[i][Jmax].nx;
		A4 = (Cell[i][Jmax - 1].T + Cell[i][Jmax - 2].T) * FaceI[i][Jmax - 1].nx / 2;
		Cell[i][Jmax - 1].T_x = (A1 - A2 + A3 - A4) / Cell[i][Jmax - 1].Vol;
		B1 = (Cell[i][Jmax - 1].T + Cell[i + 1][Jmax - 1].T) * FaceJ[i + 1][Jmax - 1].ny / 2;
		B2 = (Cell[i][Jmax - 1].T + Cell[i - 1][Jmax - 1].T) * FaceJ[i][Jmax - 1].ny / 2;
		B3 = FaceI[i][Jmax].T * FaceI[i][Jmax].ny;
		B4 = (Cell[i][Jmax - 1].T + Cell[i][Jmax - 2].T) * FaceI[i][Jmax - 1].ny / 2;
		Cell[i][Jmax - 1].T_y = (B1 - B2 + B3 - B4) / Cell[i][Jmax - 1].Vol;
				
		A1 = (Cell[i][1].u + Cell[i + 1][1].u) * FaceJ[i + 1][1].nx / 2;
		A2 = (Cell[i][1].u + Cell[i - 1][1].u) * FaceJ[i][1].nx / 2;
		A3 = (Cell[i][1].u + Cell[i][2].u) * FaceI[i][2].nx / 2;
		A4 = FaceI[i][1].u * FaceI[i][1].nx;
		Cell[i][1].u_x = (A1 - A2 + A3 - A4) / Cell[i][1].Vol;
		B1 = (Cell[i][1].u + Cell[i + 1][1].u) * FaceJ[i + 1][1].ny / 2;
		B2 = (Cell[i][1].u + Cell[i - 1][1].u) * FaceJ[i][1].ny / 2;
		B3 = (Cell[i][1].u + Cell[i][2].u) * FaceI[i][2].ny / 2;
		B4 = FaceI[i][1].u * FaceI[i][1].ny;
		Cell[i][1].u_y = (B1 - B2 + B3 - B4) / Cell[i][1].Vol;
		A1 = (Cell[i][Jmax - 1].u + Cell[i + 1][Jmax - 1].u) * FaceJ[i + 1][Jmax - 1].nx / 2;
		A2 = (Cell[i][Jmax - 1].u + Cell[i - 1][Jmax - 1].u) * FaceJ[i][Jmax - 1].nx / 2;
		A3 = FaceI[i][Jmax].u * FaceI[i][Jmax].nx;
		A4 = (Cell[i][Jmax - 1].u + Cell[i][Jmax - 2].u) * FaceI[i][Jmax - 1].nx / 2;
		Cell[i][Jmax - 1].u_x = (A1 - A2 + A3 - A4) / Cell[i][Jmax - 1].Vol;
		B1 = (Cell[i][Jmax - 1].u + Cell[i + 1][Jmax - 1].u) * FaceJ[i + 1][Jmax - 1].ny / 2;
		B2 = (Cell[i][Jmax - 1].u + Cell[i - 1][Jmax - 1].u) * FaceJ[i][Jmax - 1].ny / 2;
		B3 = FaceI[i][Jmax].u * FaceI[i][Jmax].ny;
		B4 = (Cell[i][Jmax - 1].u + Cell[i][Jmax - 2].u) * FaceI[i][Jmax - 1].ny / 2;
		Cell[i][Jmax - 1].u_y = (B1 - B2 + B3 - B4) / Cell[i][Jmax - 1].Vol;

		A1 = (Cell[i][1].v + Cell[i + 1][1].v) * FaceJ[i + 1][1].nx / 2;
		A2 = (Cell[i][1].v + Cell[i - 1][1].v) * FaceJ[i][1].nx / 2;
		A3 = (Cell[i][1].v + Cell[i][2].v) * FaceI[i][2].nx / 2;
		A4 = FaceI[i][1].v * FaceI[i][1].nx;
		Cell[i][1].v_x = (A1 - A2 + A3 - A4) / Cell[i][1].Vol;
		B1 = (Cell[i][1].v + Cell[i + 1][1].v) * FaceJ[i + 1][1].ny / 2;
		B2 = (Cell[i][1].v + Cell[i - 1][1].v) * FaceJ[i][1].ny / 2;
		B3 = (Cell[i][1].v + Cell[i][2].v) * FaceI[i][2].ny / 2;
		B4 = FaceI[i][1].v * FaceI[i][1].ny;
		Cell[i][1].v_y = (B1 - B2 + B3 - B4) / Cell[i][1].Vol;
		A1 = (Cell[i][Jmax - 1].v + Cell[i + 1][Jmax - 1].v) * FaceJ[i + 1][Jmax - 1].nx / 2;
		A2 = (Cell[i][Jmax - 1].v + Cell[i - 1][Jmax - 1].v) * FaceJ[i][Jmax - 1].nx / 2;
		A3 = FaceI[i][Jmax].v * FaceI[i][Jmax].nx;
		A4 = (Cell[i][Jmax - 1].v + Cell[i][Jmax - 2].v) * FaceI[i][Jmax - 1].nx / 2;
		Cell[i][Jmax - 1].v_x = (A1 - A2 + A3 - A4) / Cell[i][Jmax - 1].Vol;
		B1 = (Cell[i][Jmax - 1].v + Cell[i + 1][Jmax - 1].v) * FaceJ[i + 1][Jmax - 1].ny / 2;
		B2 = (Cell[i][Jmax - 1].v + Cell[i - 1][Jmax - 1].v) * FaceJ[i][Jmax - 1].ny / 2;
		B3 = FaceI[i][Jmax].v * FaceI[i][Jmax].ny;
		B4 = (Cell[i][Jmax - 1].v + Cell[i][Jmax - 2].v) * FaceI[i][Jmax - 1].ny / 2;
		Cell[i][Jmax - 1].v_y = (B1 - B2 + B3 - B4) / Cell[i][Jmax - 1].Vol;
	}
	for (j = 2; j <= Jmax - 2; j++)
	{
		FaceJ[1][j].T = 0.5 * (Cell[1][j].T + DummyJW[1][j].T);
		FaceJ[Imax][j].T = 0.5 * (Cell[Imax - 1][j].T + DummyJE[1][j].T);
		FaceJ[1][j].u = 0.5 * (Cell[1][j].u + DummyJW[1][j].u);
		FaceJ[Imax][j].u = 0.5 * (Cell[Imax - 1][j].u + DummyJE[1][j].u);
		FaceJ[1][j].v = 0.5 * (Cell[1][j].v + DummyJW[1][j].v);
		FaceJ[Imax][j].v = 0.5 * (Cell[Imax - 1][j].v + DummyJE[1][j].v);

		A1 = (Cell[1][j].T + Cell[2][j].T) * FaceJ[2][j].nx / 2;
		A2 = FaceJ[1][j].T * FaceJ[1][j].nx;
		A3 = (Cell[1][j].T + Cell[1][j + 1].T) * FaceI[1][j + 1].nx / 2;
		A4 = (Cell[1][j].T + Cell[1][j - 1].T) * FaceI[1][j].nx / 2;
		Cell[1][j].T_x = (A1 - A2 + A3 - A4) / Cell[1][j].Vol;
		B1 = (Cell[1][j].T + Cell[2][j].T) * FaceJ[2][j].ny / 2;
		B2 = FaceJ[1][j].T * FaceJ[1][j].ny;
		B3 = (Cell[1][j].T + Cell[1][j + 1].T) * FaceI[1][j + 1].ny / 2;
		B4 = (Cell[1][j].T + Cell[1][j - 1].T) * FaceI[1][j].ny / 2;
		Cell[1][j].T_y = (B1 - B2 + B3 - B4) / Cell[1][j].Vol;
		A1 = FaceJ[Imax][j].T * FaceJ[Imax][j].nx;
		A2 = (Cell[Imax - 1][j].T + Cell[Imax - 2][j].T) * FaceJ[Imax - 1][j].nx / 2;
		A3 = (Cell[Imax - 1][j].T + Cell[Imax - 1][j + 1].T) * FaceI[Imax - 1][j + 1].nx / 2;
		A4 = (Cell[Imax - 1][j].T + Cell[Imax - 1][j - 1].T) * FaceI[Imax - 1][j].nx / 2;
		Cell[Imax - 1][j].T_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][j].Vol;
		B1 = FaceJ[Imax][j].T * FaceJ[Imax][j].ny;
		B2 = (Cell[Imax - 1][j].T + Cell[Imax - 2][j].T) * FaceJ[Imax - 1][j].ny / 2;
		B3 = (Cell[Imax - 1][j].T + Cell[Imax - 1][j + 1].T) * FaceI[Imax - 1][j + 1].ny / 2;
		B4 = (Cell[Imax - 1][j].T + Cell[Imax - 1][j - 1].T) * FaceI[Imax - 1][j].ny / 2;
		Cell[Imax - 1][j].T_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][j].Vol;

		A1 = (Cell[1][j].u + Cell[2][j].u) * FaceJ[2][j].nx / 2;
		A2 = FaceJ[1][j].u * FaceJ[1][j].nx;
		A3 = (Cell[1][j].u + Cell[1][j + 1].u) * FaceI[1][j + 1].nx / 2;
		A4 = (Cell[1][j].u + Cell[1][j - 1].u) * FaceI[1][j].nx / 2;
		Cell[1][j].u_x = (A1 - A2 + A3 - A4) / Cell[1][j].Vol;
		B1 = (Cell[1][j].u + Cell[2][j].u) * FaceJ[2][j].ny / 2;
		B2 = FaceJ[1][j].u * FaceJ[1][j].ny;
		B3 = (Cell[1][j].u + Cell[1][j + 1].u) * FaceI[1][j + 1].ny / 2;
		B4 = (Cell[1][j].u + Cell[1][j - 1].u) * FaceI[1][j].ny / 2;
		Cell[1][j].u_y = (B1 - B2 + B3 - B4) / Cell[1][j].Vol;
		A1 = FaceJ[Imax][j].u * FaceJ[Imax][j].nx;
		A2 = (Cell[Imax - 1][j].u + Cell[Imax - 2][j].u) * FaceJ[Imax - 1][j].nx / 2;
		A3 = (Cell[Imax - 1][j].u + Cell[Imax - 1][j + 1].u) * FaceI[Imax - 1][j + 1].nx / 2;
		A4 = (Cell[Imax - 1][j].u + Cell[Imax - 1][j - 1].u) * FaceI[Imax - 1][j].nx / 2;
		Cell[Imax - 1][j].u_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][j].Vol;
		B1 = FaceJ[Imax][j].u * FaceJ[Imax][j].ny;
		B2 = (Cell[Imax - 1][j].u + Cell[Imax - 2][j].u) * FaceJ[Imax - 1][j].ny / 2;
		B3 = (Cell[Imax - 1][j].u + Cell[Imax - 1][j + 1].u) * FaceI[Imax - 1][j + 1].ny / 2;
		B4 = (Cell[Imax - 1][j].u + Cell[Imax - 1][j - 1].u) * FaceI[Imax - 1][j].ny / 2;
		Cell[Imax - 1][j].u_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][j].Vol;

		A1 = (Cell[1][j].v + Cell[2][j].v) * FaceJ[2][j].nx / 2;
		A2 = FaceJ[1][j].v * FaceJ[1][j].nx;
		A3 = (Cell[1][j].v + Cell[1][j + 1].v) * FaceI[1][j + 1].nx / 2;
		A4 = (Cell[1][j].v + Cell[1][j - 1].v) * FaceI[1][j].nx / 2;
		Cell[1][j].v_x = (A1 - A2 + A3 - A4) / Cell[1][j].Vol;
		B1 = (Cell[1][j].v + Cell[2][j].v) * FaceJ[2][j].ny / 2;
		B2 = FaceJ[1][j].v * FaceJ[1][j].ny;
		B3 = (Cell[1][j].v + Cell[1][j + 1].v) * FaceI[1][j + 1].ny / 2;
		B4 = (Cell[1][j].v + Cell[1][j - 1].v) * FaceI[1][j].ny / 2;
		Cell[1][j].v_y = (B1 - B2 + B3 - B4) / Cell[1][j].Vol;
		A1 = FaceJ[Imax][j].v * FaceJ[Imax][j].nx;
		A2 = (Cell[Imax - 1][j].v + Cell[Imax - 2][j].v) * FaceJ[Imax - 1][j].nx / 2;
		A3 = (Cell[Imax - 1][j].v + Cell[Imax - 1][j + 1].v) * FaceI[Imax - 1][j + 1].nx / 2;
		A4 = (Cell[Imax - 1][j].v + Cell[Imax - 1][j - 1].v) * FaceI[Imax - 1][j].nx / 2;
		Cell[Imax - 1][j].v_x = (A1 - A2 + A3 - A4) / Cell[Imax - 1][j].Vol;
		B1 = FaceJ[Imax][j].v * FaceJ[Imax][j].ny;
		B2 = (Cell[Imax - 1][j].v + Cell[Imax - 2][j].v) * FaceJ[Imax - 1][j].ny / 2;
		B3 = (Cell[Imax - 1][j].v + Cell[Imax - 1][j + 1].v) * FaceI[Imax - 1][j + 1].ny / 2;
		B4 = (Cell[Imax - 1][j].v + Cell[Imax - 1][j - 1].v) * FaceI[Imax - 1][j].ny / 2;
		Cell[Imax - 1][j].v_y = (B1 - B2 + B3 - B4) / Cell[Imax - 1][j].Vol;
	}
	
	//3.内部
	for (i = 2; i <= Imax - 2; i++)
	{
		for (j = 2; j <= Jmax - 2; j++)
		{
			A1 = (Cell[i][j].T + Cell[i + 1][j].T) * FaceJ[i + 1][j].nx / 2;
			A2 = (Cell[i][j].T + Cell[i - 1][j].T) * FaceJ[i][j].nx / 2;
			A3 = (Cell[i][j].T + Cell[i][j + 1].T) * FaceI[i][j + 1].nx / 2;
			A4 = (Cell[i][j].T + Cell[i][j - 1].T) * FaceI[i][j].nx / 2;
			Cell[i][j].T_x = (A1 - A2 + A3 - A4) / Cell[i][j].Vol;
			B1 = (Cell[i][j].T + Cell[i + 1][j].T) * FaceJ[i + 1][j].ny / 2;
			B2 = (Cell[i][j].T + Cell[i - 1][j].T) * FaceJ[i][j].ny / 2;
			B3 = (Cell[i][j].T + Cell[i][j + 1].T) * FaceI[i][j + 1].ny / 2;
			B4 = (Cell[i][j].T + Cell[i][j - 1].T) * FaceI[i][j].ny / 2;
			Cell[i][j].T_y = (B1 - B2 + B3 - B4) / Cell[i][j].Vol;

			A1 = (Cell[i][j].u + Cell[i + 1][j].u) * FaceJ[i + 1][j].nx / 2;
			A2 = (Cell[i][j].u + Cell[i - 1][j].u) * FaceJ[i][j].nx / 2;
			A3 = (Cell[i][j].u + Cell[i][j + 1].u) * FaceI[i][j + 1].nx / 2;
			A4 = (Cell[i][j].u + Cell[i][j - 1].u) * FaceI[i][j].nx / 2;
			Cell[i][j].u_x = (A1 - A2 + A3 - A4) / Cell[i][j].Vol;
			B1 = (Cell[i][j].u + Cell[i + 1][j].u) * FaceJ[i + 1][j].ny / 2;
			B2 = (Cell[i][j].u + Cell[i - 1][j].u) * FaceJ[i][j].ny / 2;
			B3 = (Cell[i][j].u + Cell[i][j + 1].u) * FaceI[i][j + 1].ny / 2;
			B4 = (Cell[i][j].u + Cell[i][j - 1].u) * FaceI[i][j].ny / 2;
			Cell[i][j].u_y = (B1 - B2 + B3 - B4) / Cell[i][j].Vol;

			A1 = (Cell[i][j].v + Cell[i + 1][j].v) * FaceJ[i + 1][j].nx / 2;
			A2 = (Cell[i][j].v + Cell[i - 1][j].v) * FaceJ[i][j].nx / 2;
			A3 = (Cell[i][j].v + Cell[i][j + 1].v) * FaceI[i][j + 1].nx / 2;
			A4 = (Cell[i][j].v + Cell[i][j - 1].v) * FaceI[i][j].nx / 2;
			Cell[i][j].v_x = (A1 - A2 + A3 - A4) / Cell[i][j].Vol;
			B1 = (Cell[i][j].v + Cell[i + 1][j].v) * FaceJ[i + 1][j].ny / 2;
			B2 = (Cell[i][j].v + Cell[i - 1][j].v) * FaceJ[i][j].ny / 2;
			B3 = (Cell[i][j].v + Cell[i][j + 1].v) * FaceI[i][j + 1].ny / 2;
			B4 = (Cell[i][j].v + Cell[i][j - 1].v) * FaceI[i][j].ny / 2;
			Cell[i][j].v_y = (B1 - B2 + B3 - B4) / Cell[i][j].Vol;
		}
	}

	//4.虚元上的偏导数，应该计算一层，总共四个边
	for (i = 1; i <= Imax - 1; i++)
	{
		DummyIS[i][1].T_x = DummyIS[i][1].T_y = 0;
		DummyIS[i][1].u_x = Cell[i][1].u_x;
		DummyIS[i][1].u_y = Cell[i][1].u_y;
		DummyIS[i][1].v_x = Cell[i][1].v_x;
		DummyIS[i][1].v_y = Cell[i][1].v_y;//物面边界
		DummyIN[i][1].T_x = DummyIN[i][1].T_y = DummyIN[i][1].u_x = DummyIN[i][1].u_y = DummyIN[i][1].v_x = DummyIN[i][1].v_y = 0;//远场边界
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		DummyJW[1][j].T_x = Cell[Imax - 1][j].T_x;
		DummyJW[1][j].T_y = Cell[Imax - 1][j].T_y;
		DummyJW[1][j].u_x = Cell[Imax - 1][j].u_x;
		DummyJW[1][j].u_y = Cell[Imax - 1][j].u_y;
		DummyJW[1][j].v_x = Cell[Imax - 1][j].v_x;
		DummyJW[1][j].v_y = Cell[Imax - 1][j].v_y;//左割缝

		DummyJE[1][j].T_x = Cell[1][j].T_x;
		DummyJE[1][j].T_y = Cell[1][j].T_y;
		DummyJE[1][j].u_x = Cell[1][j].u_x;
		DummyJE[1][j].u_y = Cell[1][j].u_y;
		DummyJE[1][j].v_x = Cell[1][j].v_x;
		DummyJE[1][j].v_y = Cell[1][j].v_y;//左割缝
	}
}
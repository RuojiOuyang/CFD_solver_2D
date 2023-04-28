/*
Convective������������ĵ�Ķ�����
*/

#include "Main.h"

using namespace std;

void Convective()
{
	//1.���б߽���غ���
	//1.1�߽紦���غ���
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][1].FU[1] = 0.5 * (Cell[i][1].U[1] + DummyIS[i][1].U[1]);
		FaceI[i][1].FU[2] = 0.5 * (Cell[i][1].U[2] + DummyIS[i][1].U[2]);
		FaceI[i][1].FU[3] = 0.5 * (Cell[i][1].U[3] + DummyIS[i][1].U[3]);
		FaceI[i][1].FU[4] = 0.5 * (Cell[i][1].U[4] + DummyIS[i][1].U[4]);

		FaceI[i][Jmax].FU[1] = 0.5 * (Cell[i][Jmax - 1].U[1] + DummyIN[i][1].U[1]);
		FaceI[i][Jmax].FU[2] = 0.5 * (Cell[i][Jmax - 1].U[2] + DummyIN[i][1].U[2]);
		FaceI[i][Jmax].FU[3] = 0.5 * (Cell[i][Jmax - 1].U[3] + DummyIN[i][1].U[3]);
		FaceI[i][Jmax].FU[4] = 0.5 * (Cell[i][Jmax - 1].U[4] + DummyIN[i][1].U[4]);
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		FaceJ[1][j].FU[1] = 0.5 * (Cell[1][j].U[1] + DummyJW[1][j].U[1]);
		FaceJ[1][j].FU[2] = 0.5 * (Cell[1][j].U[2] + DummyJW[1][j].U[2]);
		FaceJ[1][j].FU[3] = 0.5 * (Cell[1][j].U[3] + DummyJW[1][j].U[3]);
		FaceJ[1][j].FU[4] = 0.5 * (Cell[1][j].U[4] + DummyJW[1][j].U[4]);

		FaceJ[Imax][j].FU[1] = 0.5 * (Cell[Imax - 1][j].U[1] + DummyJE[1][j].U[1]);
		FaceJ[Imax][j].FU[2] = 0.5 * (Cell[Imax - 1][j].U[2] + DummyJE[1][j].U[2]);
		FaceJ[Imax][j].FU[3] = 0.5 * (Cell[Imax - 1][j].U[3] + DummyJE[1][j].U[3]);
		FaceJ[Imax][j].FU[4] = 0.5 * (Cell[Imax - 1][j].U[4] + DummyJE[1][j].U[4]);
	}

	//1.2�ڲ��߽紦���غ���
	for (j = 2; j <= Jmax - 1; j++)//i=1:Imax-1,j=2:Jmax-1��Ӧ�ı߶����ڲ�
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			FaceI[i][j].FU[1] = 0.5 * (Cell[i][j - 1].U[1] + Cell[i][j].U[1]);//�߽��ϵĵ����������ĵ���ƽ��
			FaceI[i][j].FU[2] = 0.5 * (Cell[i][j - 1].U[2] + Cell[i][j].U[2]);
			FaceI[i][j].FU[3] = 0.5 * (Cell[i][j - 1].U[3] + Cell[i][j].U[3]);
			FaceI[i][j].FU[4] = 0.5 * (Cell[i][j - 1].U[4] + Cell[i][j].U[4]);
		}
	}
	for (j = 1; j <= Jmax - 1; j++)//i=2:Imax-1,j=1:Jmax-1��Ӧ�ı߶����ڲ�
	{
		for (i = 2; i <= Imax - 1; i++)
		{
			FaceJ[i][j].FU[1] = 0.5 * (Cell[i - 1][j].U[1] + Cell[i][j].U[1]);//�߽��ϵĵ����������ĵ���ƽ��
			FaceJ[i][j].FU[2] = 0.5 * (Cell[i - 1][j].U[2] + Cell[i][j].U[2]);
			FaceJ[i][j].FU[3] = 0.5 * (Cell[i - 1][j].U[3] + Cell[i][j].U[3]);
			FaceJ[i][j].FU[4] = 0.5 * (Cell[i - 1][j].U[4] + Cell[i][j].U[4]);
		}
	}

	//2.�ɱ߽��ϵ��غ�������߽��ϵ�ͨ��
	double A, B;
	for (j = 1; j <= Jmax; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			A = FaceI[i][j].FU[1] * (pow(FaceI[i][j].FU[2] / FaceI[i][j].FU[1], 2.0) + pow(FaceI[i][j].FU[3] / FaceI[i][j].FU[1], 2.0)) / 2.0; //A���ܶ�*����
			B = (k - 1.0) * (FaceI[i][j].FU[4] - A);//B��ѹ��
			FaceI[i][j].FF[1] = FaceI[i][j].FU[1] * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny);
			FaceI[i][j].FF[2] = FaceI[i][j].FU[2] * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny) + B * FaceI[i][j].nx;
			FaceI[i][j].FF[3] = FaceI[i][j].FU[3] * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny) + B * FaceI[i][j].ny;
			FaceI[i][j].FF[4] = (FaceI[i][j].FU[4] + B) * (FaceI[i][j].FU[2] / FaceI[i][j].FU[1] * FaceI[i][j].nx + FaceI[i][j].FU[3] / FaceI[i][j].FU[1] * FaceI[i][j].ny);
		}
	}
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax; i++)
		{
			A = FaceJ[i][j].FU[1] * (pow(FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1], 2.0) + pow(FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1], 2.0)) / 2.0;
			B = (k - 1.0) * (FaceJ[i][j].FU[4] - A);
			FaceJ[i][j].FF[1] = FaceJ[i][j].FU[1] * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny);
			FaceJ[i][j].FF[2] = FaceJ[i][j].FU[2] * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny) + B * FaceJ[i][j].nx;
			FaceJ[i][j].FF[3] = FaceJ[i][j].FU[3] * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny) + B * FaceJ[i][j].ny;
			FaceJ[i][j].FF[4] = (FaceJ[i][j].FU[4] + B) * (FaceJ[i][j].FU[2] / FaceJ[i][j].FU[1] * FaceJ[i][j].nx + FaceJ[i][j].FU[3] / FaceJ[i][j].FU[1] * FaceJ[i][j].ny);

			//if (FaceJ[i][j].FF[1] > 10)
				//cout << "error" << endl;
		}
	}

	//2.�������б߽��ͨ���󣬿����嵥Ԫ��ͨ����Ҳ����Qij
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].QC[1] = FaceI[i][j + 1].FF[1] - FaceI[i][j].FF[1] + FaceJ[i + 1][j].FF[1] - FaceJ[i][j].FF[1];
			Cell[i][j].QC[2] = FaceI[i][j + 1].FF[2] - FaceI[i][j].FF[2] + FaceJ[i + 1][j].FF[2] - FaceJ[i][j].FF[2];
			Cell[i][j].QC[3] = FaceI[i][j + 1].FF[3] - FaceI[i][j].FF[3] + FaceJ[i + 1][j].FF[3] - FaceJ[i][j].FF[3];
			Cell[i][j].QC[4] = FaceI[i][j + 1].FF[4] - FaceI[i][j].FF[4] + FaceJ[i + 1][j].FF[4] - FaceJ[i][j].FF[4];
		}
	}
}
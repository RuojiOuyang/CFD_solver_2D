/*
Local_time_step�������㵱��ʱ�䲽
*/

#include "Main.h"

using namespace std;

void Min_time_step()
{
	double A, B, C, D, E, F, G, L;//���㵱��ʱ�䲽ʱ���м�������ֲ�����
	double mintime;

	mintime = 100;

	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			A = 0.5 * (FaceJ[i][j].nx + FaceJ[i + 1][j].nx);//���ʽ�Ӿ��ǰ�ʽ(2-1-13)�����һ���SI���x����ļ���һ����
			B = 0.5 * (FaceJ[i][j].ny + FaceJ[i + 1][j].ny);//���ʽ�Ӿ��ǰ�ʽ(2-1-13)�����һ���SI���y����ļ���һ����
			C = 0.5 * (FaceI[i][j].nx + FaceI[i][j + 1].nx);//���ʽ�Ӿ��ǰ�ʽ(2-1-13)����ڶ����SI���x����ļ���һ����
			D = 0.5 * (FaceI[i][j].ny + FaceI[i][j + 1].ny);//���ʽ�Ӿ��ǰ�ʽ(2-1-13)����ڶ����SI���y����ļ���һ����
			E = fabs(Cell[i][j].u * A + Cell[i][j].v * B);//ʽ(2-1-13)�����һ��
			F = fabs(Cell[i][j].u * C + Cell[i][j].v * D);//ʽ(2-1-13)����ڶ���
			G = fabs(sqrt(A * A + B * B));
			L = fabs(sqrt(C * C + D * D));//ʽ(2-1-13)���������
			Cell[i][j].localdt = CFL * Cell[i][j].Vol / (E + F + Cell[i][j].c * (G + L));//���������ʱ�䲽���������ĵ㴦

			if (Cell[i][j].localdt < mintime)//��¼�µ�ǰ��Сʱ�䲽
			{
				mintime = Cell[i][j].localdt;
			}
		}
	}
	
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].dt = mintime;//�ƽ�ʱ�������Сʱ�䲽
		}
	}

	totaltime += mintime;
}
/*
Initialization����������ʼ������
������������������ʼ��
*/

#include "Main.h"

using namespace std;

void Initialization()
{
	//������ֵ������ʼ�����������������ֵ��Ϊ����ֵ
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			Cell[i][j].Den = Dent;
			Cell[i][j].P = Pt;
			Cell[i][j].u = ut;
			Cell[i][j].v = vt;
			Cell[i][j].E = Et;
			Cell[i][j].H = Ht;
			Cell[i][j].T = Tt;
			Cell[i][j].c = ct;
			Cell[i][j].Ma = sqrt(ut * ut + vt * ct) / ct;
			Cell[i][j].miu = miut;
			Cell[i][j].miubl = miublt;
		}
	}
}
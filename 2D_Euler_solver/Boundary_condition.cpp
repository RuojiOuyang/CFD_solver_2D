/*
Boundary_condition���������߽紦��ֵ���߽��ϵ�ͨ��Face[][].FF������Convective�������Qij
*/

#include "Main.h"

using namespace std;

void Boundary_condition()
{
	//1.Զ���߽�����
	//i=1;Imax-1,j=JmaxΪԶ����Զ���Ƕ�Face���Ե�
	double c, qtinf, qninf, qti, qni, qt, qn;
	double norm_S;//����ʸ����ģ
	double qn_if, a_if;//�жϽ��������ǳ����ٵ���
	double Re_z, Rinf_f;

	for (i = 1; i <= Imax - 1; i++)
	{
		norm_S = sqrt(FaceI[i][Jmax].nx * FaceI[i][Jmax].nx + FaceI[i][Jmax].ny * FaceI[i][Jmax].ny);//����ʸ����ģ

		qtinf = (ut * FaceI[i][Jmax].ny - vt * FaceI[i][Jmax].nx) / norm_S;//�����������ٶ�
		qninf = (ut * FaceI[i][Jmax].nx + vt * FaceI[i][Jmax].ny) / norm_S;//�����ķ����ٶ�
		qti = (Cell[i][Jmax - 1].u * FaceI[i][Jmax].ny - Cell[i][Jmax - 1].v * FaceI[i][Jmax].nx) / norm_S;//�ڳ��������ٶ�
		qni = (Cell[i][Jmax - 1].u * FaceI[i][Jmax].nx + Cell[i][Jmax - 1].v * FaceI[i][Jmax].ny) / norm_S;//�ڳ��ķ����ٶ�
		
		if (sqrt(ut * ut + vt * vt) < ct)//������
		{
			Re_z = qni + 2 * Cell[i][Jmax - 1].c / (k - 1.0);
			Rinf_f = qninf - 2 * ct / (k - 1.0);

			qn_if = 1.0 / 2.0 * (Re_z + Rinf_f);
			a_if = (k - 1.0) / 4 * (Re_z - Rinf_f);

			if (qn_if <= 0)//�����ٽ���
			{
				//��ȡ����
				//�����ٶ�ȡ����
				qt = qtinf;
				//�����ٶ������ⳡ����
				qn = 1.0 / 2.0 * (Re_z + Rinf_f);
				//���������ⳡ����
				c = (k - 1.0) / 4 * (Re_z - Rinf_f);
				//���������������߽��ϵ�Den��u��v��E
				FaceI[i][Jmax].Den = pow(pow(Dent, k) * c * c / k / Pt, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
			else//�����ٳ���
			{
				//��ȡ�ڳ�
				//�����ٶ�ȡ�ڳ�
				qt = qti;
				//�����ٶ������ⳡ����
				qn = 1.0 / 2.0 * (Re_z + Rinf_f);
				//���������ⳡ����
				c = (k - 1.0) / 4 * (Re_z - Rinf_f);
				//���������������߽��ϵ�Den��u��v��E
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den, k) * c * c / k / Cell[i][Jmax - 1].P, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
		}
		else//������
		{
			if (ut * FaceI[i][Jmax].nx + vt * FaceI[i][Jmax].ny < 0)//�����ٽ���
			{
				//��ȡ����
				//�����ٶ�ȡ����
				qt = qtinf;
				//�����ٶ�����������
				qn = 1.0 / 2.0 * (qninf + 2 * ct / (k - 1.0) + qninf - 2 * ct / (k - 1.0));
				//��������������
				c = (k - 1.0) / 4 * (qninf + 2 * ct / (k - 1.0) - qninf + 2 * ct / (k - 1.0));
				//���������������߽��ϵ�Den��u��v��E
				FaceI[i][Jmax].Den = pow(pow(Dent, k) * c * c / k / Pt, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
			else//�����ٳ���
			{
				//��ȡ�ڳ�
				//�����ٶ�ȡ�ڳ�
				qt = qti;
				//�����ٶ����ڳ�����
				qn = 1.0 / 2.0 * (qni + 2 * Cell[i][Jmax - 1].c / (k - 1.0) + qni - 2 * Cell[i][Jmax - 1].c / (k - 1.0));
				//�������ڳ�����
				c = (k - 1) / 4 * (qni + 2 * Cell[i][Jmax - 1].c / (k - 1) - qni + 2 * Cell[i][Jmax - 1].c / (k - 1));
				//���������������߽��ϵ�Den��u��v��E
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den, k) * c * c / k / Cell[i][Jmax - 1].P, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
		}
	}
	//1.2Զ���߽�����غ���
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][Jmax].FU[1] = FaceI[i][Jmax].Den;
		FaceI[i][Jmax].FU[2] = FaceI[i][Jmax].Den * FaceI[i][Jmax].u;
		FaceI[i][Jmax].FU[3] = FaceI[i][Jmax].Den * FaceI[i][Jmax].v;
		FaceI[i][Jmax].FU[4] = FaceI[i][Jmax].Den * FaceI[i][Jmax].E;
	}
	//1.3��Զ���߽���غ�������ͨ��
	double A, B;
	for (i = 1; i <= Imax - 1; i++)
	{
		A = FaceI[i][Jmax].FU[1] * (pow(FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1], 2.0) + pow(FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1], 2.0)) / 2.0;//A��den*����
		B = (k - 1.0) * (FaceI[i][Jmax].FU[4] - A);//B��ѹ��
		FaceI[i][Jmax].FF[1] = FaceI[i][Jmax].FU[1] * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny);
		FaceI[i][Jmax].FF[2] = FaceI[i][Jmax].FU[2] * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny) + B * FaceI[i][Jmax].nx;
		FaceI[i][Jmax].FF[3] = FaceI[i][Jmax].FU[3] * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny) + B * FaceI[i][Jmax].ny;
		FaceI[i][Jmax].FF[4] = (FaceI[i][Jmax].FU[4] + B) * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny);
	}

	//2.����߽�����
	//i=1:Imax-1,j=1Ϊ����
	double S1, S2, S3;
	//2.1��ֵ�õ�����߽磨�±߽磩����
	for (i = 1; i <= Imax - 1; i++)
	{
		S3 = sqrt(pow((Cell[i][3].x - 0.5 * Node[i][1].x - 0.5 * Node[i + 1][1].x), 2) + pow((Cell[i][3].y - 0.5 * Node[i][1].y - 0.5 * Node[i + 1][1].y), 2));
		S2 = sqrt(pow((Cell[i][2].x - 0.5 * Node[i][1].x - 0.5 * Node[i + 1][1].x), 2) + pow((Cell[i][2].y - 0.5 * Node[i][1].y - 0.5 * Node[i + 1][1].y), 2));
		S1 = sqrt(pow((Cell[i][1].x - 0.5 * Node[i][1].x - 0.5 * Node[i + 1][1].x), 2) + pow((Cell[i][1].y - 0.5 * Node[i][1].y - 0.5 * Node[i + 1][1].y), 2));
		FaceI[i][1].P = (S3 / (S3 - S1) * S2 / (S2 - S1)) * Cell[i][1].P - (S3 / (S3 - S2) * S1 / (S2 - S1)) * Cell[i][2].P + (S2 / (S3 - S2) * S1 / (S3 - S1)) * Cell[i][3].P;
		FaceI[i][1].Den = (S3 / (S3 - S1) * S2 / (S2 - S1)) * Cell[i][1].Den - (S3 / (S3 - S2) * S1 / (S2 - S1)) * Cell[i][2].Den + (S2 / (S3 - S2) * S1 / (S3 - S1)) * Cell[i][3].Den;
		FaceI[i][1].T = (S3 / (S3 - S1) * S2 / (S2 - S1)) * Cell[i][1].T - (S3 / (S3 - S2) * S1 / (S2 - S1)) * Cell[i][2].T + (S2 / (S3 - S2) * S1 / (S3 - S1)) * Cell[i][3].T;
		FaceI[i][1].E = (S3 / (S3 - S1) * S2 / (S2 - S1)) * Cell[i][1].E - (S3 / (S3 - S2) * S1 / (S2 - S1)) * Cell[i][2].E + (S2 / (S3 - S2) * S1 / (S3 - S1)) * Cell[i][3].E;
		FaceI[i][1].u = sqrt(2 * (FaceI[i][1].E - FaceI[i][1].P / FaceI[i][1].Den / (k - 1.0)));
		FaceI[i][1].v = 0;
		//cout << sqrt(pow((Cell[i][Jmax - 1].x - 0.5 * Node[i][1].x - 0.5 * Node[i + 1][1].x), 2) + pow((Cell[i][Jmax - 1].y - 0.5 * Node[i][1].y - 0.5 * Node[i + 1][1].y), 2)) << endl;
	}
	//2.2�������߽磨�±߽磩ͨ��
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][1].FF[1] = 0;
		FaceI[i][1].FF[2] = FaceI[i][1].P * FaceI[i][1].nx;
		FaceI[i][1].FF[3] = FaceI[i][1].P * FaceI[i][1].ny;
		FaceI[i][1].FF[4] = 0;
	}

	//3.���߽�����
	//3.1�õ����ұ߽���غ���
	for (j = 1; j <= Jmax - 1; j++)//i=Imax����i=1��һ��
	{
		FaceJ[Imax][j].FU[1] = 1.0 / 2.0 * (Cell[Imax - 1][j].U[1] + Cell[1][j].U[1]);
		FaceJ[Imax][j].FU[2] = 1.0 / 2.0 * (Cell[Imax - 1][j].U[2] + Cell[1][j].U[2]);
		FaceJ[Imax][j].FU[3] = 1.0 / 2.0 * (Cell[Imax - 1][j].U[3] + Cell[1][j].U[3]);
		FaceJ[Imax][j].FU[4] = 1.0 / 2.0 * (Cell[Imax - 1][j].U[4] + Cell[1][j].U[4]);
		FaceJ[1][j].FU[1] = FaceJ[Imax][j].FU[1];
		FaceJ[1][j].FU[2] = FaceJ[Imax][j].FU[2];
		FaceJ[1][j].FU[3] = FaceJ[Imax][j].FU[3];
		FaceJ[1][j].FU[4] = FaceJ[Imax][j].FU[4];
	}
	//3.2�ɸ��߽���غ�������ͨ��
	for (j = 1; j <= Jmax - 1; j++)
	{
		A = FaceJ[Imax][j].FU[1] * (pow(FaceJ[Imax][j].FU[2] / FaceJ[Imax][j].FU[1], 2.0) + pow(FaceJ[Imax][j].FU[3] / FaceJ[Imax][j].FU[1], 2.0)) / 2.0;//A��den*����
		B = (k - 1.0) * (FaceJ[Imax][j].FU[4] - A);//B��ѹ��
		FaceJ[Imax][j].FF[1] = FaceJ[Imax][j].FU[1] * (FaceJ[Imax][j].FU[2] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].nx + FaceJ[Imax][j].FU[3] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].ny);
		FaceJ[Imax][j].FF[2] = FaceJ[Imax][j].FU[2] * (FaceJ[Imax][j].FU[2] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].nx + FaceJ[Imax][j].FU[3] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].ny) + B * FaceJ[Imax][j].nx;
		FaceJ[Imax][j].FF[3] = FaceJ[Imax][j].FU[3] * (FaceJ[Imax][j].FU[2] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].nx + FaceJ[Imax][j].FU[3] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].ny) + B * FaceJ[Imax][j].ny;
		FaceJ[Imax][j].FF[4] = FaceJ[Imax][j].FU[4] * (FaceJ[Imax][j].FU[2] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].nx + FaceJ[Imax][j].FU[3] / FaceJ[Imax][j].FU[1] * FaceJ[Imax][j].ny);

		A = FaceJ[1][j].FU[1] * (pow(FaceJ[1][j].FU[2] / FaceJ[1][j].FU[1], 2.0) + pow(FaceJ[1][j].FU[3] / FaceJ[1][j].FU[1], 2.0)) / 2.0;//A��den*����
		B = (k - 1.0) * (FaceJ[1][j].FU[4] - A);//B��ѹ��
		FaceJ[1][j].FF[1] = FaceJ[1][j].FU[1] * (FaceJ[1][j].FU[2] / FaceJ[1][j].FU[1] * FaceJ[1][j].nx + FaceJ[1][j].FU[3] / FaceJ[1][j].FU[1] * FaceJ[1][j].ny);
		FaceJ[1][j].FF[2] = FaceJ[1][j].FU[2] * (FaceJ[1][j].FU[2] / FaceJ[1][j].FU[1] * FaceJ[1][j].nx + FaceJ[1][j].FU[3] / FaceJ[1][j].FU[1] * FaceJ[1][j].ny) + B * FaceJ[1][j].nx;
		FaceJ[1][j].FF[3] = FaceJ[1][j].FU[3] * (FaceJ[1][j].FU[2] / FaceJ[1][j].FU[1] * FaceJ[1][j].nx + FaceJ[1][j].FU[3] / FaceJ[1][j].FU[1] * FaceJ[1][j].ny) + B * FaceJ[1][j].ny;
		FaceJ[1][j].FF[4] = FaceJ[1][j].FU[4] * (FaceJ[1][j].FU[2] / FaceJ[1][j].FU[1] * FaceJ[1][j].nx + FaceJ[1][j].FU[3] / FaceJ[1][j].FU[1] * FaceJ[1][j].ny);
	}
}
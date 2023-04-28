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

		qtinf = (utinf * FaceI[i][Jmax].ny - vtinf * FaceI[i][Jmax].nx) / norm_S;//�����������ٶ�
		qninf = (utinf * FaceI[i][Jmax].nx + vtinf * FaceI[i][Jmax].ny) / norm_S;//�����ķ����ٶ�
		qti = (Cell[i][Jmax - 1].u * ctinf * FaceI[i][Jmax].ny - Cell[i][Jmax - 1].v * ctinf * FaceI[i][Jmax].nx) / norm_S;//�ڳ��������ٶ�
		qni = (Cell[i][Jmax - 1].u * ctinf * FaceI[i][Jmax].nx + Cell[i][Jmax - 1].v * ctinf * FaceI[i][Jmax].ny) / norm_S;//�ڳ��ķ����ٶ�

		if (sqrt(ut * ut + vt * vt) < ct)//������
		{
			Re_z = qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1.0);
			Rinf_f = qninf - 2 * ct * ctinf / (k - 1.0);

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
				FaceI[i][Jmax].Den = pow(pow(Dentinf, k) * c * c / k / Ptinf, 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R / Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf/ ctinf;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
				//cout << FaceI[i][Jmax].E << endl;
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
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den * Dentinf, k) * c * c / k / (Cell[i][Jmax - 1].P * Dentinf * ctinf * ctinf), 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R / Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf / ctinf;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
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
				qn = 1.0 / 2.0 * (qninf + 2 * ctinf / (k - 1.0) + qninf - 2 * ctinf / (k - 1.0));
				//��������������
				c = (k - 1.0) / 4 * (qninf + 2 * ctinf / (k - 1.0) - qninf + 2 * ctinf / (k - 1.0));
				//���������������߽��ϵ�Den��u��v��E
				FaceI[i][Jmax].Den = pow(pow(Dentinf, k) * c * c / k / Ptinf, 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R /Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf / ctinf;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
			else//�����ٳ���
			{
				//��ȡ�ڳ�
				//�����ٶ�ȡ�ڳ�
				qt = qti;
				//�����ٶ����ڳ�����
				qn = 1.0 / 2.0 * (qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1.0) + qni - 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1.0));
				//�������ڳ�����
				c = (k - 1) / 4 * (qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1) - qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1));
				//���������������߽��ϵ�Den��u��v��E
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den * Dentinf, k) * c * c / k / (Cell[i][Jmax - 1].P * Dentinf * ctinf * ctinf), 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R / Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf / ctinf;;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
		}
	}
	/*
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
	*/
}
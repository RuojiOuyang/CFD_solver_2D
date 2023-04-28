/*
BoundaryCondition���������߽��ϵ�ͨ��������Ϊʱ���
*/

#include "Main.h"

using namespace std;
using namespace arma;

void BoundaryCondition()
{
	//1.Զ���߽�����
	//i=1;Imax-1,j=JmaxΪԶ����Զ���Ƕ�Face���Ե�
	mat norm_S;//����ʸ����ģ
	mat qn_inf, Rinf_f;
	mat qni, Re_z;
	mat qn_b, c_b;
	mat Zerosmat;//ȫ0������
	mat Idenmat;//ȫ1������
	umat ZZ1, ZZ2;
	mat Den1, P1, T1, u1, v1, E1;
	mat Den2, P2, T2, u2, v2, E2;

	/*mat qt_inf, qti;
	double qt, qn, c;

	Idenmat = ones((long long)Imax - 1, 1);

	norm_S = arma::sqrt(FaceI_nx.col((long long)Jmax - 1) % FaceI_nx.col((long long)Jmax - 1) + FaceI_ny.col((long long)Jmax - 1) % FaceI_ny.col((long long)Jmax - 1));//����ʸ����ģ

	qt_inf = (ut * FaceI_ny.col((long long)Jmax - 1) - vt * FaceI_nx.col((long long)Jmax - 1)) / norm_S;//�����������ٶ�
	qn_inf = (ut * FaceI_nx.col((long long)Jmax - 1) + vt * FaceI_ny.col((long long)Jmax - 1)) / norm_S;//�����ķ����ٶ�

	Rinf_f = qn_inf - 2 * ct * Idenmat / (k - 1.0);

	qti = (Cell_u.col((long long)Jmax - 2) % FaceI_ny.col((long long)Jmax - 1) - Cell_v.col((long long)Jmax - 2) % FaceI_nx.col((long long)Jmax - 1)) / norm_S;//�ڳ��������ٶ�
	qni = (Cell_u.col((long long)Jmax - 2) % FaceI_nx.col((long long)Jmax - 1) + Cell_v.col((long long)Jmax - 2) % FaceI_ny.col((long long)Jmax - 1)) / norm_S;//�ڳ��ķ����ٶ�
	Re_z = qni + 2 * Cell_c.col((long long)Jmax - 2) / (k - 1.0);
	qn_b = 1.0 / 2.0 * (Re_z + Rinf_f);

	for (i = 0; i < Imax - 1; i++)
	{
		if (qn_b(i, 0) <= 0)//�����ٽ���
		{
			qt = qt_inf(i, 0);//�����ٶ�ȡ����
			qn = 1.0 / 2.0 * (Re_z(i, 0) + Rinf_f(i, 0));//�����ٶ������ⳡ����
			c = (k - 1.0) / 4 * (Re_z(i, 0) - Rinf_f(i, 0));//���������ⳡ����
			//���������������߽��ϵ�Den��u��v��E
			FaceI_Den(i, (long long)Jmax - 1) = pow(pow(Dent, k) * c * c / k / Pt, 1.0 / (k - 1.0));
			FaceI_T(i, (long long)Jmax - 1) = c * c / k / R;
			FaceI_P(i, (long long)Jmax - 1) = FaceI_Den(i, (long long)Jmax - 1) * R * FaceI_T(i, (long long)Jmax - 1);
			FaceI_u(i, (long long)Jmax - 1) = (qt * FaceI_ny(i, (long long)Jmax - 1) + qn * FaceI_nx(i, (long long)Jmax - 1)) / norm_S(i, 0);
			FaceI_v(i, (long long)Jmax - 1) = (-qt * FaceI_nx(i, (long long)Jmax - 1) + qn * FaceI_ny(i, (long long)Jmax - 1)) / norm_S(i, 0);
			FaceI_E(i, (long long)Jmax - 1) = FaceI_P(i, (long long)Jmax - 1) / (k - 1.0) / FaceI_Den(i, (long long)Jmax - 1) +
				(FaceI_u(i, (long long)Jmax - 1) * FaceI_u(i, (long long)Jmax - 1) + FaceI_v(i, (long long)Jmax - 1) * FaceI_v(i, (long long)Jmax - 1)) / 2.0;
		}
		else//�����ٳ���
		{
			qt = qti(i, 0);//�����ٶ�ȡ�ڳ�
			qn = 1.0 / 2.0 * (Re_z(i, 0) + Rinf_f(i, 0));//�����ٶ������ⳡ����
			c = (k - 1.0) / 4 * (Re_z(i, 0) - Rinf_f(i, 0));//���������ⳡ����
			//���������������߽��ϵ�Den��u��v��E
			FaceI_Den(i, (long long)Jmax - 1) = pow(pow(Cell_Den(i, (long long)Jmax - 2), k) * c * c / k / Cell_P(i, (long long)Jmax - 2), 1.0 / (k - 1.0));
			FaceI_T(i, (long long)Jmax - 1) = c * c / k / R;
			FaceI_P(i, (long long)Jmax - 1) = FaceI_Den(i, (long long)Jmax - 1) * R * FaceI_T(i, (long long)Jmax - 1);
			FaceI_u(i, (long long)Jmax - 1) = (qt * FaceI_ny(i, (long long)Jmax - 1) + qn * FaceI_nx(i, (long long)Jmax - 1)) / norm_S(i, 0);
			FaceI_v(i, (long long)Jmax - 1) = (-qt * FaceI_nx(i, (long long)Jmax - 1) + qn * FaceI_ny(i, (long long)Jmax - 1)) / norm_S(i, 0);
			FaceI_E(i, (long long)Jmax - 1) = FaceI_P(i, (long long)Jmax - 1) / (k - 1.0) / FaceI_Den(i, (long long)Jmax - 1) +
				(FaceI_u(i, (long long)Jmax - 1) * FaceI_u(i, (long long)Jmax - 1) + FaceI_v(i, (long long)Jmax - 1) * FaceI_v(i, (long long)Jmax - 1)) / 2.0;
		}
	}*/

	Zerosmat = zeros((long long)Imax - 1, 1);
	Idenmat = ones((long long)Imax - 1, 1);

	norm_S = arma::sqrt(FaceI_nx.col((long long)Jmax - 1) % FaceI_nx.col((long long)Jmax - 1) + FaceI_ny.col((long long)Jmax - 1) % FaceI_ny.col((long long)Jmax - 1));//����ʸ����ģ

	qn_inf = (ut * FaceI_nx.col((long long)Jmax - 1) + vt * FaceI_ny.col((long long)Jmax - 1)) / norm_S;//�����ķ����ٶ�
	Rinf_f = qn_inf - 2 * ct / (k - 1.0);//�ⳡ������������

	qni = (Cell_u.col((long long)Jmax - 2) % FaceI_nx.col((long long)Jmax - 1) + Cell_v.col((long long)Jmax - 2) % FaceI_ny.col((long long)Jmax - 1)) / norm_S;//�ڳ��ķ����ٶ�
	Re_z = qni + 2 * Cell_c.col((long long)Jmax - 2) / (k - 1.0);//�ڳ�������������

	qn_b = 1.0 / 2.0 * (Re_z + Rinf_f);//�߽�ķ����ٶ�
	c_b = (k - 1.0) / 4.0 * (Re_z - Rinf_f);//�߽������

	ZZ1 = (qn_b < Zerosmat);
	ZZ2 = (qn_b >= Zerosmat);

	Den1 = Dent * arma::pow(c_b % c_b / ct / ct, 1.0 / (k - 1.0));
	P1 = Den1 % c_b % c_b / k;
	T1 = P1 / Den1 / R;
	u1 = ut + FaceI_nx.col((long long)Jmax - 1) % (qn_b - qn_inf) / norm_S;
	v1 = vt + FaceI_ny.col((long long)Jmax - 1) % (qn_b - qn_inf) / norm_S;
	E1 = P1 / (k - 1.0) / Den1 + (u1 % u1 + v1 % v1) / 2.0;

	Den2 = Cell_Den.col((long long)Jmax - 2) % arma::pow(c_b % c_b / Cell_c.col((long long)Jmax - 2) / Cell_c.col((long long)Jmax - 2), 1.0 / (k - 1.0));
	P2 = Den2 % c_b % c_b / k;
	T2 = P2 / Den2 / R;
	u2 = Cell_u.col((long long)Jmax - 2) + FaceI_nx.col((long long)Jmax - 1) % (qn_b - qn_inf) / norm_S;
	v2 = Cell_v.col((long long)Jmax - 2) + FaceI_ny.col((long long)Jmax - 1) % (qn_b - qn_inf) / norm_S;
	E2 = P2 / (k - 1.0) / Den2 + (u2 % u2 + v2 % v2) / 2.0;

	FaceI_Den.col((long long)Jmax - 1) = Den1 % ZZ1 + Den2 % ZZ2;
	FaceI_P.col((long long)Jmax - 1) = P1 % ZZ1 + P2 % ZZ2;
	FaceI_T.col((long long)Jmax - 1) = T1 % ZZ1 + T2 % ZZ2;
	FaceI_u.col((long long)Jmax - 1) = u1 % ZZ1 + u2 % ZZ2;
	FaceI_v.col((long long)Jmax - 1) = v1 % ZZ1 + v2 % ZZ2;
	FaceI_E.col((long long)Jmax - 1) = E1 % ZZ1 + E2 % ZZ2;

	//1.2Զ���߽�����غ���
	FaceI_U1.col((long long)Jmax - 1) = FaceI_Den.col((long long)Jmax - 1);
	FaceI_U2.col((long long)Jmax - 1) = FaceI_Den.col((long long)Jmax - 1) % FaceI_u.col((long long)Jmax - 1);
	FaceI_U3.col((long long)Jmax - 1) = FaceI_Den.col((long long)Jmax - 1) % FaceI_v.col((long long)Jmax - 1);
	FaceI_U4.col((long long)Jmax - 1) = FaceI_Den.col((long long)Jmax - 1) % FaceI_E.col((long long)Jmax - 1);
	//1.3��Զ���߽���غ�������ͨ��
	mat AA, BB;
	AA = FaceI_U1.col((long long)Jmax - 1) % (arma::pow(FaceI_U2.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1), 2.0) +
		arma::pow(FaceI_U3.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1), 2.0)) / 2.0;//A��den*����
	BB = (k - 1.0) * (FaceI_U4.col((long long)Jmax - 1) - AA);//B��ѹ��
	FaceI_F1.col((long long)Jmax - 1) = FaceI_U1.col((long long)Jmax - 1) % (FaceI_U2.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) %
		FaceI_nx.col((long long)Jmax - 1) + FaceI_U3.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) % FaceI_ny.col((long long)Jmax - 1));
	FaceI_F2.col((long long)Jmax - 1) = FaceI_U2.col((long long)Jmax - 1) % (FaceI_U2.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) %
		FaceI_nx.col((long long)Jmax - 1) + FaceI_U3.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) % FaceI_ny.col((long long)Jmax - 1)) +
		BB % FaceI_nx.col((long long)Jmax - 1);
	FaceI_F3.col((long long)Jmax - 1) = FaceI_U3.col((long long)Jmax - 1) % (FaceI_U2.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) %
		FaceI_nx.col((long long)Jmax - 1) + FaceI_U3.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) % FaceI_ny.col((long long)Jmax - 1)) +
		BB % FaceI_ny.col((long long)Jmax - 1);
	FaceI_F4.col((long long)Jmax - 1) = (FaceI_U4.col((long long)Jmax - 1) + BB) % (FaceI_U2.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) %
		FaceI_nx.col((long long)Jmax - 1) + FaceI_U3.col((long long)Jmax - 1) / FaceI_U1.col((long long)Jmax - 1) % FaceI_ny.col((long long)Jmax - 1));

	//2.����߽�����
	//i=1:Imax-1,j=1Ϊ����
	mat S1, S2, S3;
	//2.1��ֵ�õ�����߽磨�±߽磩����
	S3 = arma::sqrt(arma::pow(Cell_x.col(2) - 0.5 * Node_x(span(0, (long long)Imax - 2), span(0)) - 0.5 * Node_x(span(1, (long long)Imax - 1), span(0)), 2.0) +
		arma::pow(Cell_y.col(2) - 0.5 * Node_y(span(0, (long long)Imax - 2), span(0)) - 0.5 * Node_y(span(1, (long long)Imax - 1), span(0)), 2));
	S2 = arma::sqrt(arma::pow(Cell_x.col(1) - 0.5 * Node_x(span(0, (long long)Imax - 2), span(0)) - 0.5 * Node_x(span(1, (long long)Imax - 1), span(0)), 2.0) +
		arma::pow(Cell_y.col(1) - 0.5 * Node_y(span(0, (long long)Imax - 2), span(0)) - 0.5 * Node_y(span(1, (long long)Imax - 1), span(0)), 2));
	S1 = arma::sqrt(arma::pow(Cell_x.col(0) - 0.5 * Node_x(span(0, (long long)Imax - 2), span(0)) - 0.5 * Node_x(span(1, (long long)Imax - 1), span(0)), 2.0) +
		arma::pow(Cell_y.col(0) - 0.5 * Node_y(span(0, (long long)Imax - 2), span(0)) - 0.5 * Node_y(span(1, (long long)Imax - 1), span(0)), 2));
	FaceI_P.col(0) = (S3 / (S3 - S1) % S2 / (S2 - S1)) % Cell_P.col(0) - (S3 / (S3 - S2) % S1 / (S2 - S1)) % Cell_P.col(1) + (S2 / (S3 - S2) % S1 / (S3 - S1)) % Cell_P.col(2);
	//2.2�������߽磨�±߽磩ͨ��
	FaceI_F1.col(0).zeros();
	FaceI_F2.col(0) = FaceI_P.col(0) % FaceI_nx.col(0);
	FaceI_F3.col(0) = FaceI_P.col(0) % FaceI_ny.col(0);
	FaceI_F4.col(0).zeros();

	//3.���߽�����
	//3.1�õ����ұ߽���غ�����Ϊ�˼���ʱ��ɱ����͵�2����ѭ��д����һ��
	FaceJ_U1.row((long long)Imax - 1) = 1.0 / 2.0 * (Cell_U1.row((long long)Imax - 2) + Cell_U1.row(0));
	FaceJ_U2.row((long long)Imax - 1) = 1.0 / 2.0 * (Cell_U2.row((long long)Imax - 2) + Cell_U2.row(0));
	FaceJ_U3.row((long long)Imax - 1) = 1.0 / 2.0 * (Cell_U3.row((long long)Imax - 2) + Cell_U3.row(0));
	FaceJ_U4.row((long long)Imax - 1) = 1.0 / 2.0 * (Cell_U4.row((long long)Imax - 2) + Cell_U4.row(0));
	FaceJ_U1.row(0) = FaceJ_U1.row((long long)Imax - 1);
	FaceJ_U2.row(0) = FaceJ_U2.row((long long)Imax - 1);
	FaceJ_U3.row(0) = FaceJ_U3.row((long long)Imax - 1);
	FaceJ_U4.row(0) = FaceJ_U4.row((long long)Imax - 1);
	//3.2�ɸ��߽���غ�������ͨ��
	AA = FaceJ_U1.row((long long)Imax - 1) % (arma::pow(FaceJ_U2.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1), 2.0) +
		arma::pow(FaceJ_U3.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1), 2.0)) / 2.0;//A��den*����
	BB = (k - 1.0) * (FaceJ_U4.row((long long)Imax - 1) - AA);//B��ѹ��
	FaceJ_F1.row((long long)Imax - 1) = FaceJ_U1.row((long long)Imax - 1) % (FaceJ_U2.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) %
		FaceJ_nx.row((long long)Imax - 1) + FaceJ_U3.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) % FaceJ_ny.row((long long)Imax - 1));
	FaceJ_F2.row((long long)Imax - 1) = FaceJ_U2.row((long long)Imax - 1) % (FaceJ_U2.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) %
		FaceJ_nx.row((long long)Imax - 1) + FaceJ_U3.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) % FaceJ_ny.row((long long)Imax - 1)) + BB % FaceJ_nx.row((long long)Imax - 1);
	FaceJ_F3.row((long long)Imax - 1) = FaceJ_U3.row((long long)Imax - 1) % (FaceJ_U2.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) %
		FaceJ_nx.row((long long)Imax - 1) + FaceJ_U3.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) % FaceJ_ny.row((long long)Imax - 1)) + BB % FaceJ_ny.row((long long)Imax - 1);
	FaceJ_F4.row((long long)Imax - 1) = (FaceJ_U4.row((long long)Imax - 1) + BB) % (FaceJ_U2.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) %
		FaceJ_nx.row((long long)Imax - 1) + FaceJ_U3.row((long long)Imax - 1) / FaceJ_U1.row((long long)Imax - 1) % FaceJ_ny.row((long long)Imax - 1));
	AA = FaceJ_U1.row(0) % (arma::pow(FaceJ_U2.row(0) / FaceJ_U1.row(0), 2.0) +
		arma::pow(FaceJ_U3.row(0) / FaceJ_U1.row(0), 2.0)) / 2.0;//A��den*����
	BB = (k - 1.0) * (FaceJ_U4.row(0) - AA);//B��ѹ��
	FaceJ_F1.row(0) = FaceJ_U1.row(0) % (FaceJ_U2.row(0) / FaceJ_U1.row(0) %
		FaceJ_nx.row(0) + FaceJ_U3.row(0) / FaceJ_U1.row(0) % FaceJ_ny.row(0));
	FaceJ_F2.row(0) = FaceJ_U2.row(0) % (FaceJ_U2.row(0) / FaceJ_U1.row(0) %
		FaceJ_nx.row(0) + FaceJ_U3.row(0) / FaceJ_U1.row(0) % FaceJ_ny.row(0)) +
		BB % FaceJ_nx.row(0);
	FaceJ_F3.row(0) = FaceJ_U3.row(0) % (FaceJ_U2.row(0) / FaceJ_U1.row(0) %
		FaceJ_nx.row(0) + FaceJ_U3.row(0) / FaceJ_U1.row(0) % FaceJ_ny.row(0)) +
		BB % FaceJ_ny.row(0);
	FaceJ_F4.row(0) = (FaceJ_U4.row(0) + BB) % (FaceJ_U2.row(0) / FaceJ_U1.row(0) %
		FaceJ_nx.row(0) + FaceJ_U3.row(0) / FaceJ_U1.row(0) % FaceJ_ny.row(0));
}
#include "Main.h"
#include "Read_mesh.h"
#include "Geometry.h"
#include "Initialization.h"
#include "Solve_conservation.h"
#include "Runge_Kutta.h"
#include "Residual.h"
#include "Result_deal.h"
#include "Putout.h"

using namespace std;

int main(int argc, char* argv[])
{
	totaltime = 0;

	ifstream read;
	read.open("yuanzhudata.txt");
	read >> Imax >> Jmax;//��ǰ��������ǰ�������ֱ���i,j��������������
	read.close();
	n_dummy = 3;//��Ԫ������
	Dynamic_array();//�����ڴ�

	//��ʼ��������
	//�����趨
	k = 1.4;//���ȱ�
	R = 287.06;//���峣��
	T0 = 288.16;//�¶Ȳο�ֵ
	Ts = 110;//Sutherland����
	P0 = 101325.0;//ѹ���ο�ֵ
	mu0 = 1.7894e-5;//��Ӧ��T0�Ķ���ճ��
	Pr = 0.71;//��������
	CFL = 2;//CFL
	AOA = 0;//���ǣ���λ��
	RKalpha[1] = 1.0 / 4.0;//RK����ϵ��
	RKalpha[2] = 1.0 / 6.0;
	RKalpha[3] = 3.0 / 8.0;
	RKalpha[4] = 0.5;
	RKalpha[5] = 1.0;
	Iternum_max = 10000;//����������

	//��������
	//L = 50 * 1e-3;//�����ߴ�
	Re = 228;//������ŵ��
	Ma = 0.0002;
	Ttinf = T0 / (1.0 + (k - 1.0) / 2.0 * Ma * Ma);//�����¶�
	Ttinf = T0;
	mutinf = mu0 * pow(Ttinf / T0, 1.5) * (T0 + Ts) / (Ttinf + Ts);//����������
	//P0 = Re * R * Ttinf * mutinf / (pow(Ttinf / T0, k / (k - 1.0)) * sqrt(k * R * Ttinf) * Ma * L);//���ѹ���ο�ֵ
	Ptinf = P0 * pow(Ttinf / T0, k / (k - 1.0));//����ѹ��
	ctinf = sqrt(k * R * Ttinf);//��������
	Dentinf = Ptinf / (R * Ttinf);//�����ܶ�
	utinf = ctinf * Ma * cos(AOA / 180 * PI);//����x�����ٶ�
	vtinf = ctinf * Ma * sin(AOA / 180 * PI);//����y�����ٶ�
	Etinf = Ptinf / (k - 1.0) / Dentinf + (utinf * utinf + vtinf * vtinf) / 2.0;//����������
	Htinf = Etinf + Ptinf / Dentinf;//��������
	L = Re * mutinf / Dentinf / Ma / ctinf;//�����ߴ�

	Tt = 1;
	Pt = Ptinf / (Dentinf * ctinf * ctinf);
	ct = 1;
	Dent = 1;
	ut = utinf / ctinf;
	vt = vtinf / ctinf;
	Et = Etinf / ctinf / ctinf;
	Ht = Htinf / ctinf / ctinf;
	mut = 1;

	//��������
	Read_mesh();//������
	Geometry();//�������񣬵õ�����������Լ����ʸ��

	//���㲿��
	Initialization();//��������ĳ�ʼ������ʼ������

	Solve_conservation();//��������ĵ��ϵ��غ���

	ofstream outx;
	outx.open("Res.dat");//д���
	for (n = 1; n <= Iternum_max; n++)//����ѭ������
	{
		Runge_Kutta();//RK����ǰ���
		Residual();//����в�
		cout << n << "   " << Resid << endl;
		outx << n << "   " << Resid << endl;
		if (Resid <= 1e-18)
		{
			break;
		}
	}
	outx.close();

	Result_deal();//��������

	Putout();//�������

	cout << totaltime * L / ctinf << endl;

	for (n = Iternum_max + 1;; n++)
	{
		Runge_Kutta();//RK����ǰ���
		Residual();//����в�
		cout << n << "   " << Resid << endl;
		outx << n << "   " << Resid << endl;
		if (n % 10000 == 0)
		{
			Result_deal();//��������
			Putout();//�������
			cout << totaltime * L / ctinf << endl;
		}
	}

	return 0;
}
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
	T0 = 288.16;//��ƽ���¶Ȳο�ֵ
	Ts = 110;//Sutherland����
	mu0 = 1.7894e-5;//��Ӧ��T0�Ķ���ճ��
	P0 = 101325.0;//��ƽ��ѹ���ο�ֵ
	c0 = 340.28;//��ƽ�����ٲο�ֵ
	cv = R / (k - 1.0);
	cp = k * cv;
	Pr = 0.71;//��������
	CFL = 2;//CFL
	AOA = 0;//���ǣ���λ��
	RKalpha[1] = 1.0 / 4.0;//RK����ϵ��
	RKalpha[2] = 1.0 / 6.0;
	RKalpha[3] = 3.0 / 8.0;
	RKalpha[4] = 0.5;
	RKalpha[5] = 1.0;
	Iternum_max = 30000;//����������

	//SA����ģ���еĳ���
	Cb1 = 0.1355;
	Cb2 = 0.622;
	delta = 2.0 / 3.0;
	writek = 0.41;
	rlim = 10;
	Cw1 = Cb1 / writek / writek + (1 + Cb2) / delta;
	Cw2 = 0.3;
	Cw3 = 2.0;
	Cv1 = 7.1;
	Prt = 0.9;//����Prantl��

	//��������
	Ma = 0.2;//���������
	Tt = T0 / (1.0 + (k - 1.0) / 2.0 * Ma * Ma);//�����¶�
	Pt = P0 * pow(Tt / T0, k / (k - 1.0));//����ѹ��
	ct = sqrt(k * R * Tt);//��������
	Dent = Pt / (R * Tt);//�����ܶ�
	ut = ct * Ma * cos(AOA / 180 * PI);//����x�����ٶ�
	vt = ct * Ma * sin(AOA / 180 * PI);//����y�����ٶ�
	Et = Pt / (k - 1.0) / Dent + (ut * ut + vt * vt) / 2.0;//����������
	Ht = Et + Pt / Dent;//��������
	miut = mu0 * pow(Tt / T0, 1.5) * (T0 + Ts) / (Tt + Ts);//��������ճ��
	miublt = 0.1 * miut / Dent;//����miuһ���ˣ�SAģ���е�

	cout << Dent * Ma * ct * 50 * 1e-7 / miut << endl;

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
		if (Resid <= 1e-28)
		{
			break;
		}
	}

	Result_deal();//��������

	Putout();//�������

	cout << totaltime << endl;

	for (n = Iternum_max + 1;; n++)
	{
		Runge_Kutta();//RK����ǰ���
		Residual();//����в�
		cout << n << "   " << Resid << endl;
		outx << n << "   " << Resid << endl;
		if (n % 1000 == 0)
		{
			Result_deal();//��������
			Putout();//�������
			cout << totaltime << endl;
		}
	}

	outx.close();
	return 0;
}
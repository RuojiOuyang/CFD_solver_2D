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
	n_dummy = 3;//��Ԫ������

	//��������
	ifstream read;
	read.open("naca0012data.txt");
	read >> Imax >> Jmax;//��ǰ��������ǰ�������ֱ���i,j��������������
	read.close();
	Dynamic_array();//�����ڴ�

	//��ʼ��������
	//�����趨
	k = 1.4;//���ȱ�
	R = 287.06;//���峣��
	T0 = 288.16;//��ƽ���¶Ȳο�ֵ
	P0 = 101325.0;//��ƽ��ѹ���ο�ֵ
	c0 = 340.28;//��ƽ�����ٲο�ֵ
	CFL = 2;//CFL
	AOA = 1.25;//���ǣ���λ��
	RKalpha[1] = 1.0 / 4.0;//RK����ϵ��
	RKalpha[2] = 1.0 / 6.0;
	RKalpha[3] = 3.0 / 8.0;
	RKalpha[4] = 0.5;
	RKalpha[5] = 1.0;
	Iternum_max = 2000;//����������

	//��������
	Ma = 0.8;//���������
	Tt = T0 / (1.0 + (k - 1.0) / 2.0 * Ma * Ma);//�����¶�
	Pt = P0 * pow(Tt / T0, k / (k - 1.0));//����ѹ��
	ct = sqrt(k * R * Tt);//��������
	Dent = Pt / (R * Tt);//�����ܶ�
	ut = ct * Ma * cos(AOA / 180 * PI);//����x�����ٶ�
	vt = ct * Ma * sin(AOA / 180 * PI);//����y�����ٶ�
	Et = Pt / (k - 1.0) / Dent + (ut * ut + vt * vt) / 2.0;//����������
	Ht = Et + Pt / Dent;//��������
	
	//��������
	Read_mesh();//������
	Geometry();//�������񣬵õ�����������Լ����ʸ��

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
		if (Resid <= 1e-10)
		{
			break;
		}
	}
	outx.close();

	Result_deal();//��������

	Putout();//�������
	
	return 0;

}
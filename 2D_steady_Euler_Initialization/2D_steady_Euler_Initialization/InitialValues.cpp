/*
InitialValues��������������ʼ״̬
*/

#include "Main.h"
#include "DynamicArray.h"
#include "ReadMesh.h"
#include "Geometry.h"

using namespace std;
using namespace arma;

void InitialValues()
{
	//1.�йؿ����ĳ���
	k = 1.4;//���ȱ�
	R = 287.06;//���峣��
	T0 = 288.16;//��ƽ���¶Ȳο�ֵ
	P0 = 101325.0;//��ƽ��ѹ���ο�ֵ
	c0 = 340.28;//��ƽ�����ٲο�ֵ

	//2.��������
	double ti;
	double Omega;
	Ma = 0.6;//���������
	Tt = T0 / (1.0 + (k - 1.0) / 2.0 * Ma * Ma);//�����¶�
	Pt = P0 * pow(Tt / T0, k / (k - 1.0));//����ѹ��
	ct = sqrt(k * R * Tt);//��������
	Dent = Pt / (R * Tt);//�����ܶ�

	//Omega = 2 * ct * Ma * 0.0808 / 1.0;
	//ti = 2 * PI / Omega * ((long long)numAOA - 1.0) / M;//ʱ�����
	//AOA = 2.89 + 2.41 * sin(Omega * ti);//���ǣ���λ��
	AOA = 3;

	ut = ct * Ma * cos(AOA / 180 * PI);//����x�����ٶ�
	vt = ct * Ma * sin(AOA / 180 * PI);//����y�����ٶ�
	Et = Pt / (k - 1.0) / Dent + (ut * ut + vt * vt) / 2.0;//����������

	//3.��������Ҫ�����ڴ棬�����ڴ�֮ǰ��ҪһЩ�ؼ�������Imax,Jmax,n_dummy
	ifstream read;
	read.open("naca0012_09_17data.txt");
	read >> Imax >> Jmax;//��ǰ��������ǰ�������ֱ���i,j��������������
	read.close();
	n_dummy = 3;
	DynamicArray();

	//���������Ҫ����
	CFL = 2;//CFL��
	RKalpha[1] = 1.0 / 4.0;//RK����ϵ��
	RKalpha[2] = 1.0 / 6.0;
	RKalpha[3] = 3.0 / 8.0;
	RKalpha[4] = 0.5;
	RKalpha[5] = 1.0;
	Iternum_max = 10000;//����������

	//��������
	ReadMesh();//������
	Geometry();//�������񣬵õ�����������Լ����ʸ��
}
/*
����������������������
��Ķ������ͬ��ͷ�ļ���
Dynamic_array������Ԥ�����ڴ�
*/

#include "Variables_declear.h"

using namespace std;

CNode** Node;//�������Ϣ
CCell** Cell;//���ĵ�
CCell** DummyIS;//��Ԫ��I����׶ˣ�south
CCell** DummyIN;//��Ԫ��I���򶥶ˣ�north
CCell** DummyJW;//��Ԫ��J������ˣ�west
CCell** DummyJE;//��Ԫ��J�����Ҷˣ�east
CFace** FaceI;//I��߽磬ָ˳��i����ı߽�
CFace** FaceJ;//J��߽磬ָ˳��j����ı߽�

int i, j;//�����ѭ������
int Imax, Jmax;//���������i������j�Ƿ���
int n_dummy;//��Ԫ��ghost cell��������
int n;//��ʾ����ѭ����

//��������ĳ���
double k;//���ȱ�
double R;//���峣��
double T0;//��ƽ���¶Ȳο�ֵ
double Ts;//Sutherland����
double mu0;//��Ӧ��T0�Ķ���ճ��
double P0;//��ƽ��ѹ���ο�ֵ
double c0;//��ƽ�����ٲο�ֵ
double cv;
double cp;
double Pr;//��������

//������̵ĳ���
double CFL;//CFL
double AOA;//���ǣ���λ��
int Iternum_max;//����������
double* RKalpha;//Runge_Kutta��ϵ��
double totaltime;//�����ƽ�����ʱ��

//��������еı���
double** Den_former;//Den_former[]���ڱ�����һ���ĸ����ܶ�,���ڲв����
double** muI;//I����ļ�����������
double** muJ;//J����ļ�����������
double Resid;//�в�

//��������
double Re;//��ŵ��
double Ttinf, Ptinf, Dentinf, utinf, vtinf, Ma, ctinf, Etinf, Htinf, mutinf;//��ʼ��������
double Tt, Pt, Dent, ut, vt, ct, Et, Ht, mut;//��ʼ��������
double L;//�����ߴ�

void Dynamic_array()
{
	//����д��ά��������*�У�Ϊ���������i������j���򱣳�һ��
	Node = new CNode * [Imax + 1];//Node��Imax*Jmaxά�ġ�+1��Ϊ���ñ�ſ��Դ�1��ʼ������һ��
	Cell = new CCell * [Imax];//Cell��(Imax-1)*(Jmax-1)ά�ġ����ĵ��������������1
	DummyIS = new CCell * [Imax];//IS��(Imax-1)*n_dummyά��
	DummyIN = new CCell * [Imax];//IN��(Imax-1)*n_dummyά��
	DummyJW = new CCell * [n_dummy + 1];//JW��n_dummy*(Jmax-1)ά��
	DummyJE = new CCell * [n_dummy + 1];//JE��n_dummy*(Jmax-1)ά��
	FaceI = new CFace * [Imax];//FaceI��(Imax-1)*Jmaxά�ıߵ�������������������һ�£�ʵ���ϣ���������4*4����I����ı���3*4�ģ�J������4*3�ģ��ܹ�24����
	FaceJ = new CFace * [Imax + 1];//FaceJ��Imax*(Jmax-1)ά�ġ��ߵ�������������������һ��
	Den_former = new double* [Imax];//Den_former��Cellά��һ�£����ڱ�����һ���ĸ����ܶ�,���ڲв����
	muI = new double* [Imax];//�����������ӵ������ֱ�����ĵ������+2*��Ԫ����һ�£�������ļ�ΪI���򣬺�����ļ�ΪJ����muI��(Imax-1)*(Jmax+2*n_dummy)ά
	muJ = new double* [Imax + 2 * n_dummy + 1];//muJ��(Imax+2*n_dummy)*(Jmax-1)ά��

	for (i = 0; i <= Imax; i++)
	{
		Node[i] = new CNode[Jmax + 1];
		FaceJ[i] = new CFace[Jmax];
	}

	for (i = 0; i <= Imax - 1; i++)
	{
		Cell[i] = new CCell[Jmax];
		FaceI[i] = new CFace[Jmax + 1];
		DummyIS[i] = new CCell[n_dummy + 1];
		DummyIN[i] = new CCell[n_dummy + 1];
		Den_former[i] = new double[Jmax];
		muI[i] = new double[Jmax + 2 * n_dummy + 1];
	}

	for (i = 0; i <= Imax + 2 * n_dummy; i++)
	{
		muJ[i] = new double[Jmax];
	}

	for (i = 0; i <= n_dummy; i++)
	{
		DummyJW[i] = new CCell[Jmax];
		DummyJE[i] = new CCell[Jmax];
	}

	RKalpha = new double[6];//RK��ϵ��
}
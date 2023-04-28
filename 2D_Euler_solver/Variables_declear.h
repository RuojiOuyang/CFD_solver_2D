#pragma once

#define PI 3.1415926535

class CCell//CCell���ʾ���ĵ�
{
public:
	double x;//���ĵ�x����
	double y;//���ĵ�y����
	double Den;//�ܶ�
	double u;//x���ٶ�
	double v;//y���ٶ�
	double E;//������
	double H;//����
	double P;//ѹ��
	double T;//�¶�
	double Ma;//�����
	double Vol;//���������
	double dt;//����ʱ�䲽��
	double c;//��������
	double U[5];//�����غ������ȡ5��Ϊ�˿��Դ�1��ʼ������һ��
	double QC[5];//������Ķ���ͨ������
	double D[5];//��������˹�ճ��ͨ������
	double Med[5];//RK����������еĲ�����ʽ2-1-12�е�R��
	double U_former[5];//���ڱ���Runge-Kutta����һ�ε�����ֵ
};

class CFace//CFace���ʾ���������
{
public:
	double nx;//�ߵķ��������x�������㵱��ʱ�䲽��ʱ���õ�
	double ny;//�ߵķ��������y����
	double Den;//�ߵ��ܶ�
	double u;//x���ٶ�
	double v;//y���ٶ�
	double E;//�ߵ�������
	double H;//�ߵ�����
	double P;//�ߵ�ѹ��
	double T;//�ߵ��¶�
	double FU[5];//FU�Ǹ��߶���ͨ����������Ȼ��������������˵�����ǰ�����ƽ���غ������ټ���ͨ���ķ�ʽ������QC�ģ��������ȼ���ͨ����ƽ��
	double FF[5];//���ߵĶ���ͨ�������ʸ���ĵ��F.*S
	double FD[5];//���ߵ��˹�ճ��ͨ��
	double lamta;//ճ�������ϵ���������˹�ճ�Ե�ʱ����
	double epslion[3];//ճ����������Ӧϵ��
	double dU[5];//�غ�����һ�ײ��
	double dddU[5];//�غ��������ײ��
};

class CNode//CNode����װ������õ�
{
public:
	double x;//�����x����
	double y;//�����y����
	double Den;//����㴦�ܶ�
	double u;//����㴦x���ٶ�
	double v;//����㴦y���ٶ�
	double E;//����㴦������
	double H;//����㴦����
	double P;//����㴦ѹ��
	double T;//����㴦�¶�
	double Ma;//����㴦�����
};

extern CNode** Node;//�������Ϣ
extern CCell** Cell;//���ĵ�
extern CCell** DummyIS;//��Ԫ��I����׶ˣ�south
extern CCell** DummyIN;//��Ԫ��I���򶥶ˣ�north
extern CCell** DummyJW;//��Ԫ��J������ˣ�west
extern CCell** DummyJE;//��Ԫ��J�����Ҷˣ�east
extern CFace** FaceI;//I��߽磬ָ˳��i����ı߽�
extern CFace** FaceJ;//J��߽磬ָ˳��j����ı߽�

extern int i, j;//�����ѭ������
extern int Imax, Jmax;//���������i������j�Ƿ���
extern int n_dummy;//��Ԫ��ghost cell��������
extern int n;//��ʾ����ѭ����
extern double k;//���ȱ�
extern double R;//���峣��
extern double T0;//��ƽ���¶Ȳο�ֵ
extern double P0;//��ƽ��ѹ���ο�ֵ
extern double c0;//��ƽ�����ٲο�ֵ
extern double CFL;//CFL
extern double AOA;//���ǣ���λ��
extern int Iternum_max;//����������
extern double* RKalpha;//Runge_Kutta��ϵ��
extern double** Den_former;
extern double** muI;//I����߽�ļ�����������
extern double** muJ;//J����߽�ļ�����������
extern double Resid;//�в�
extern double Tt, Pt, Dent, ut, vt, Ma, ct, Et, Ht;//��ʼ��������

extern void Dynamic_array();
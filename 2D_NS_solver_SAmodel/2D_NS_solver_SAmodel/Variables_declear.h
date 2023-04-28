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
	double localdt;//����ʱ�䲽��
	double dt;//��Сʱ�䲽��
	double c;//��������
	double U[6];//�����غ������ȡ6��Ϊ�˿��Դ�1��ʼ������һ��
	double QC[6];//������Ķ���ͨ�����֣�����ճͨ�����֣�
	double QV[6];//���������ɢͨ�����֣���ճ��ͨ�����֣�
	double D[6];//��������˹�ճ��ͨ������
	double Med[6];//RK����������еĲ�����ʽ2-1-12�е�R��
	double U_former[6];//���ڱ���Runge-Kutta����һ�ε�����ֵ
	double T_x;//T��x�����ƫ����
	double T_y;//T��y�����ƫ����
	double u_x;//�ٶ�u��x�����ƫ����
	double u_y;//�ٶ�u��y�����ƫ����
	double v_x;//�ٶ�v��x�����ƫ����
	double v_y;//�ٶ�v��y�����ƫ����
	double Touxx;//������Touxx��
	double Touxy;//������Touxy��
	double Touyy;//������Touyy��
	double Touhx;//ʽ��3-1-1���е�Touhx
	double Touhy;//ʽ��3-1-1���е�Touhy
	double miu;//���嶯��ճ��
	double miut;//��ճճ�ȣ�����fv1֮��ģ���Den*miubl*fv1
	double miubl;//����ճ�ȣ�SAģ���еĳ�����Сv��
	double miubl_x;//����ճ����x�����ƫ����
	double miubl_y;//����ճ����y�����ƫ����
	double Touvx;//����ģ��Touvx��
	double Touvy;//����ģ��Touvy��
	double S;//Դ��
	double writeX;//��д��X��������ճ�������ճ��֮��
	double fv1;//SAģ��ϵ��
	double d;//���ĵ�����������ľ��룬SAģ��
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
	double FU[6];//FU�Ǹ��߶���ͨ����������Ȼ��������������˵�����ǰ�����ƽ���غ������ټ���ͨ���ķ�ʽ������QC�ģ��������ȼ���ͨ����ƽ��
	double FF[6];//���ߵĶ���ͨ�������ʸ���ĵ��F.*S
	double FK[6];//���ߵ���ɢͨ��
	double FD[6];//���ߵ��˹�ճ��ͨ��
	double lamta;//ճ�������ϵ���������˹�ճ�Ե�ʱ����
	double epslion[3];//ճ����������Ӧϵ��
	double dU[6];//�غ�����һ�ײ��
	double dddU[6];//�غ��������ײ��
	double miubl;//����ճ�ȣ�SAģ���е�
	double miubl_x;//����ճ�ȵ�ƫ������SAģ���е�
	double miubl_y;
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
	double v_x;
	double u_y;//��������������������
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
extern double Ts;//Sutherland����
extern double mu0;//��Ӧ��T0�Ķ���ճ��
extern double P0;//��ƽ��ѹ���ο�ֵ
extern double c0;//��ƽ�����ٲο�ֵ
extern double cv;
extern double cp;
extern double Pr;//��������
extern double CFL;//CFL
extern double AOA;//���ǣ���λ��
extern int Iternum_max;//����������
extern double* RKalpha;//Runge_Kutta��ϵ��
extern double totaltime;//�����ƽ�����ʱ��
extern double** Den_former;
extern double** muI;//I����߽�ļ�����������
extern double** muJ;//J����߽�ļ�����������
extern double Resid;//�в�
extern double Tt, Pt, Dent, ut, vt, Ma, ct, Et, Ht, miut, miublt;//��ʼ��������
extern double Cb1, Cb2, delta, writek, rlim, Cw1, Cw2, Cw3, Cv1, Prt;//SA����ģ���еĳ���

extern void Dynamic_array();
#pragma once

#include "Main.h"

#define PI 3.1415926535

//1.���ĵ�
extern arma::mat Cell_Den;//�ܶ�
extern arma::mat Cell_u;//x���ٶ�
extern arma::mat Cell_v;//y���ٶ�
extern arma::mat Cell_E;//������
extern arma::mat Cell_P;//ѹ��
extern arma::mat Cell_T;//�¶�
extern arma::mat Cell_Ma;//�����
extern arma::mat Cell_dt;//ʱ�䲽��
extern arma::mat Cell_c;//��������
extern arma::mat Cell_U1;//ͨ����һ��
extern arma::mat Cell_U2;//ͨ���ڶ���
extern arma::mat Cell_U3;//ͨ��������
extern arma::mat Cell_U4;//ͨ��������
extern arma::mat Cell_U1_former;//ǰһ��ͨ����һ��
extern arma::mat Cell_U2_former;//ǰһ��ͨ���ڶ���
extern arma::mat Cell_U3_former;//ǰһ��ͨ��������
extern arma::mat Cell_U4_former;//ǰһ��ͨ��������
extern arma::mat Cell_QC1;//������Ķ���ͨ�����ֵ�һ��
extern arma::mat Cell_QC2;//������Ķ���ͨ�����ֵڶ���
extern arma::mat Cell_QC3;//������Ķ���ͨ�����ֵ�����
extern arma::mat Cell_QC4;//������Ķ���ͨ�����ֵ�����
extern arma::mat Cell_D1;//��������˹�ճ��ͨ�����ֵ�һ��
extern arma::mat Cell_D2;//��������˹�ճ��ͨ�����ֵڶ���
extern arma::mat Cell_D3;//��������˹�ճ��ͨ�����ֵ�����
extern arma::mat Cell_D4;//��������˹�ճ��ͨ�����ֵ�����
extern arma::mat Cell_Med1;//RK����������еĲ�����һ��
extern arma::mat Cell_Med2;//RK����������еĲ����ڶ���
extern arma::mat Cell_Med3;//RK����������еĲ���������
extern arma::mat Cell_Med4;//RK����������еĲ���������
extern arma::mat Cell_x;//���ĵ�x����
extern arma::mat Cell_y;//���ĵ�y����
extern arma::mat Cell_Vol;//���������

//2.���������
//2.1I����棨���ŵģ�
extern arma::mat FaceI_Den;//�ܶ�
extern arma::mat FaceI_u;//x���ٶ�
extern arma::mat FaceI_v;//y���ٶ�
extern arma::mat FaceI_E;//������
extern arma::mat FaceI_P;//ѹ��
extern arma::mat FaceI_T;//�¶�
extern arma::mat FaceI_U1;//�����غ�����һ��
extern arma::mat FaceI_U2;//�����غ����ڶ���
extern arma::mat FaceI_U3;//�����غ���������
extern arma::mat FaceI_U4;//�����غ���������
extern arma::mat FaceI_F1;//���ߵĶ���ͨ�������ʸ���ĵ��F*S��һ��
extern arma::mat FaceI_F2;//���ߵĶ���ͨ�������ʸ���ĵ��F*S�ڶ���
extern arma::mat FaceI_F3;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
extern arma::mat FaceI_F4;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
extern arma::mat FaceI_lamta;////ճ�������ϵ���������˹�ճ�Ե�ʱ����
extern arma::mat FaceI_epslion1;//ճ����������Ӧϵ��
extern arma::mat FaceI_epslion2;//ճ����������Ӧϵ��
extern arma::mat FaceI_dU1;//�غ�����һ�ײ�ֵ�һ��
extern arma::mat FaceI_dU2;//�غ�����һ�ײ�ֵڶ���
extern arma::mat FaceI_dU3;//�غ�����һ�ײ�ֵ�����
extern arma::mat FaceI_dU4;//�غ�����һ�ײ�ֵ�����
extern arma::mat FaceI_dddU1;//�غ��������ײ�ֵ�һ��
extern arma::mat FaceI_dddU2;//�غ��������ײ�ֵڶ���
extern arma::mat FaceI_dddU3;//�غ��������ײ�ֵ�����
extern arma::mat FaceI_dddU4;//�غ��������ײ�ֵ�����
extern arma::mat FaceI_D1;//���ߵ��˹�ճ��ͨ����һ��
extern arma::mat FaceI_D2;//���ߵ��˹�ճ��ͨ���ڶ���
extern arma::mat FaceI_D3;//���ߵ��˹�ճ��ͨ��������
extern arma::mat FaceI_D4;//���ߵ��˹�ճ��ͨ��������
extern arma::mat FaceI_nx;//�ߵķ��������x����
extern arma::mat FaceI_ny;//�ߵķ��������y����
//2.2J����棨���ŵģ�
extern arma::mat FaceJ_Den;//�ܶ�
extern arma::mat FaceJ_u;//x���ٶ�
extern arma::mat FaceJ_v;//y���ٶ�
extern arma::mat FaceJ_E;//������
extern arma::mat FaceJ_P;//ѹ��
extern arma::mat FaceJ_T;//�¶�
extern arma::mat FaceJ_U1;//�����غ�����һ��
extern arma::mat FaceJ_U2;//�����غ����ڶ���
extern arma::mat FaceJ_U3;//�����غ���������
extern arma::mat FaceJ_U4;//�����غ���������
extern arma::mat FaceJ_F1;//���ߵĶ���ͨ�������ʸ���ĵ��F*S��һ��
extern arma::mat FaceJ_F2;//���ߵĶ���ͨ�������ʸ���ĵ��F*S�ڶ���
extern arma::mat FaceJ_F3;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
extern arma::mat FaceJ_F4;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
extern arma::mat FaceJ_lamta;////ճ�������ϵ���������˹�ճ�Ե�ʱ����
extern arma::mat FaceJ_epslion1;//ճ����������Ӧϵ��
extern arma::mat FaceJ_epslion2;//ճ����������Ӧϵ��
extern arma::mat FaceJ_dU1;//�غ�����һ�ײ�ֵ�һ��
extern arma::mat FaceJ_dU2;//�غ�����һ�ײ�ֵڶ���
extern arma::mat FaceJ_dU3;//�غ�����һ�ײ�ֵ�����
extern arma::mat FaceJ_dU4;//�غ�����һ�ײ�ֵ�����
extern arma::mat FaceJ_dddU1;//�غ��������ײ�ֵ�һ��
extern arma::mat FaceJ_dddU2;//�غ��������ײ�ֵڶ���
extern arma::mat FaceJ_dddU3;//�غ��������ײ�ֵ�����
extern arma::mat FaceJ_dddU4;//�غ��������ײ�ֵ�����
extern arma::mat FaceJ_D1;//���ߵ��˹�ճ��ͨ����һ��
extern arma::mat FaceJ_D2;//���ߵ��˹�ճ��ͨ���ڶ���
extern arma::mat FaceJ_D3;//���ߵ��˹�ճ��ͨ��������
extern arma::mat FaceJ_D4;//���ߵ��˹�ճ��ͨ��������
extern arma::mat FaceJ_nx;//�ߵķ��������x����
extern arma::mat FaceJ_ny;//�ߵķ��������x����

//3.�����
extern arma::mat Node_Den;//�ܶ�
extern arma::mat Node_u;//x���ٶ�
extern arma::mat Node_v;//y���ٶ�
extern arma::mat Node_E;//������
extern arma::mat Node_P;//ѹ��
extern arma::mat Node_T;//�¶�
extern arma::mat Node_Ma;//�����
extern arma::mat Node_x;
extern arma::mat Node_y;

//4.��Ԫ
//4.1I����׶ˣ�south
extern arma::mat DummyIS_Den;//�ܶ�
extern arma::mat DummyIS_u;//x���ٶ�
extern arma::mat DummyIS_v;//y���ٶ�
extern arma::mat DummyIS_E;//������
extern arma::mat DummyIS_P;//ѹ��
extern arma::mat DummyIS_T;//�¶�
extern arma::mat DummyIS_Ma;//�����
extern arma::mat DummyIS_dt;//ʱ�䲽��
extern arma::mat DummyIS_c;//��������
extern arma::mat DummyIS_U1;//ͨ����һ��
extern arma::mat DummyIS_U2;//ͨ���ڶ���
extern arma::mat DummyIS_U3;//ͨ��������
extern arma::mat DummyIS_U4;//ͨ��������
//4.2I���򶥶ˣ�north
extern arma::mat DummyIN_Den;//�ܶ�
extern arma::mat DummyIN_u;//x���ٶ�
extern arma::mat DummyIN_v;//y���ٶ�
extern arma::mat DummyIN_E;//������
extern arma::mat DummyIN_P;//ѹ��
extern arma::mat DummyIN_T;//�¶�
extern arma::mat DummyIN_Ma;//�����
extern arma::mat DummyIN_dt;//ʱ�䲽��
extern arma::mat DummyIN_c;//��������
extern arma::mat DummyIN_U1;//ͨ����һ��
extern arma::mat DummyIN_U2;//ͨ���ڶ���
extern arma::mat DummyIN_U3;//ͨ��������
extern arma::mat DummyIN_U4;//ͨ��������
//4.3J������ˣ�west
extern arma::mat DummyJW_Den;//�ܶ�
extern arma::mat DummyJW_u;//x���ٶ�
extern arma::mat DummyJW_v;//y���ٶ�
extern arma::mat DummyJW_E;//������
extern arma::mat DummyJW_P;//ѹ��
extern arma::mat DummyJW_T;//�¶�
extern arma::mat DummyJW_Ma;//�����
extern arma::mat DummyJW_dt;//ʱ�䲽��
extern arma::mat DummyJW_c;//��������
extern arma::mat DummyJW_U1;//ͨ����һ��
extern arma::mat DummyJW_U2;//ͨ���ڶ���
extern arma::mat DummyJW_U3;//ͨ��������
extern arma::mat DummyJW_U4;//ͨ��������
//4.4J�����Ҷˣ�east
extern arma::mat DummyJE_Den;//�ܶ�
extern arma::mat DummyJE_u;//x���ٶ�
extern arma::mat DummyJE_v;//y���ٶ�
extern arma::mat DummyJE_E;//������
extern arma::mat DummyJE_P;//ѹ��
extern arma::mat DummyJE_T;//�¶�
extern arma::mat DummyJE_Ma;//�����
extern arma::mat DummyJE_dt;//ʱ�䲽��
extern arma::mat DummyJE_c;//��������
extern arma::mat DummyJE_U1;//ͨ����һ��
extern arma::mat DummyJE_U2;//ͨ���ڶ���
extern arma::mat DummyJE_U3;//ͨ��������
extern arma::mat DummyJE_U4;//ͨ��������

//������Ϣ����
extern int Imax, Jmax;//���������i������j�Ƿ���

//������Ϣ����
extern int n_dummy;//��Ԫ��ghost cell��������

//�����������Ҫ�ı���
extern int i, j;//�����ѭ������
extern double CFL;//CFL��
extern int n;//��ǰ������
extern int Iternum_max;//����������
extern double* RKalpha;//Runge_Kutta��ϵ��
extern double Resid;//�в�
extern arma::mat P_former;//Den_former[]���ڱ�����һ���ĸ����ܶ�,���ڲв����
extern arma::mat muI;//I����ļ����������ӣ��˹�ճ����
extern arma::mat muJ;//J����ļ����������ӣ��˹�ճ����
extern int numAOA;//���Ǽ���
extern int M;//�����

//��������ĳ���
extern double k;//���ȱ�
extern double R;//���峣��
extern double T0;//��ƽ���¶Ȳο�ֵ
extern double P0;//��ƽ��ѹ���ο�ֵ
extern double c0;//��ƽ�����ٲο�ֵ

//��������
extern double Tt;//�����¶�
extern double Pt;//����ѹ��
extern double Dent;//�����ܶ�
extern double ut;//����ˮƽ�ٶ�
extern double vt;//������ֱ�ٶ�
extern double Ma;//���������
extern double ct;//��������
extern double Et;//��������
extern double AOA;//�������ǣ���λ��
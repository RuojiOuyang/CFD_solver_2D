/*
����������������������
��Ķ������ͬ��ͷ�ļ���
*/

/*
����������������������
��Ķ������ͬ��ͷ�ļ���
*/

#include "Main.h"
#include "VariablesDeclear.h"

//1.���ĵ�
arma::mat Cell_Den;//�ܶ�
arma::mat Cell_u;//x���ٶ�
arma::mat Cell_v;//y���ٶ�
arma::mat Cell_E;//������
arma::mat Cell_P;//ѹ��
arma::mat Cell_T;//�¶�
arma::mat Cell_Ma;//�����
arma::mat Cell_dt;//ʱ�䲽��
arma::mat Cell_c;//��������
arma::mat Cell_U1;//�غ�����һ��
arma::mat Cell_U2;//�غ����ڶ���
arma::mat Cell_U3;//�غ���������
arma::mat Cell_U4;//�غ���������
arma::mat Cell_U1_former;//ǰһ���غ�����һ��
arma::mat Cell_U2_former;//ǰһ���غ����ڶ���
arma::mat Cell_U3_former;//ǰһ���غ���������
arma::mat Cell_U4_former;//ǰһ���غ���������
arma::mat Cell_QC1;//������Ķ���ͨ�����ֵ�һ��
arma::mat Cell_QC2;//������Ķ���ͨ�����ֵڶ���
arma::mat Cell_QC3;//������Ķ���ͨ�����ֵ�����
arma::mat Cell_QC4;//������Ķ���ͨ�����ֵ�����
arma::mat Cell_D1;//��������˹�ճ��ͨ�����ֵ�һ��
arma::mat Cell_D2;//��������˹�ճ��ͨ�����ֵڶ���
arma::mat Cell_D3;//��������˹�ճ��ͨ�����ֵ�����
arma::mat Cell_D4;//��������˹�ճ��ͨ�����ֵ�����
arma::mat Cell_Med1;//RK����������еĲ�����һ��
arma::mat Cell_Med2;//RK����������еĲ����ڶ���
arma::mat Cell_Med3;//RK����������еĲ���������
arma::mat Cell_Med4;//RK����������еĲ���������
arma::mat Cell_x;//���ĵ�x����
arma::mat Cell_y;//���ĵ�y����
arma::mat Cell_Vol;//���������

//2.���������
//2.1I����棨���ŵģ�
arma::mat FaceI_Den;//�ܶ�
arma::mat FaceI_u;//x���ٶ�
arma::mat FaceI_v;//y���ٶ�
arma::mat FaceI_E;//������
arma::mat FaceI_P;//ѹ��
arma::mat FaceI_T;//�¶�
arma::mat FaceI_U1;//�����غ�����һ��
arma::mat FaceI_U2;//�����غ����ڶ���
arma::mat FaceI_U3;//�����غ���������
arma::mat FaceI_U4;//�����غ���������
arma::mat FaceI_F1;//���ߵĶ���ͨ�������ʸ���ĵ��F*S��һ��
arma::mat FaceI_F2;//���ߵĶ���ͨ�������ʸ���ĵ��F*S�ڶ���
arma::mat FaceI_F3;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
arma::mat FaceI_F4;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
arma::mat FaceI_lamta;////ճ�������ϵ���������˹�ճ�Ե�ʱ����
arma::mat FaceI_epslion1;//ճ����������Ӧϵ��
arma::mat FaceI_epslion2;//ճ����������Ӧϵ��
arma::mat FaceI_dU1;//�غ�����һ�ײ�ֵ�һ��
arma::mat FaceI_dU2;//�غ�����һ�ײ�ֵڶ���
arma::mat FaceI_dU3;//�غ�����һ�ײ�ֵ�����
arma::mat FaceI_dU4;//�غ�����һ�ײ�ֵ�����
arma::mat FaceI_dddU1;//�غ��������ײ�ֵ�һ��
arma::mat FaceI_dddU2;//�غ��������ײ�ֵڶ���
arma::mat FaceI_dddU3;//�غ��������ײ�ֵ�����
arma::mat FaceI_dddU4;//�غ��������ײ�ֵ�����
arma::mat FaceI_D1;//���ߵ��˹�ճ��ͨ����һ��
arma::mat FaceI_D2;//���ߵ��˹�ճ��ͨ���ڶ���
arma::mat FaceI_D3;//���ߵ��˹�ճ��ͨ��������
arma::mat FaceI_D4;//���ߵ��˹�ճ��ͨ��������
arma::mat FaceI_nx;//�ߵķ��������x����
arma::mat FaceI_ny;//�ߵķ��������y����
//2.2J����棨���ŵģ�
arma::mat FaceJ_Den;//�ܶ�
arma::mat FaceJ_u;//x���ٶ�
arma::mat FaceJ_v;//y���ٶ�
arma::mat FaceJ_E;//������
arma::mat FaceJ_P;//ѹ��
arma::mat FaceJ_T;//�¶�
arma::mat FaceJ_U1;//�����غ�����һ��
arma::mat FaceJ_U2;//�����غ����ڶ���
arma::mat FaceJ_U3;//�����غ���������
arma::mat FaceJ_U4;//�����غ���������
arma::mat FaceJ_F1;//���ߵĶ���ͨ�������ʸ���ĵ��F*S��һ��
arma::mat FaceJ_F2;//���ߵĶ���ͨ�������ʸ���ĵ��F*S�ڶ���
arma::mat FaceJ_F3;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
arma::mat FaceJ_F4;//���ߵĶ���ͨ�������ʸ���ĵ��F*S������
arma::mat FaceJ_lamta;//ճ�������ϵ���������˹�ճ�Ե�ʱ����
arma::mat FaceJ_epslion1;//ճ����������Ӧϵ��
arma::mat FaceJ_epslion2;//ճ����������Ӧϵ��
arma::mat FaceJ_dU1;//�غ�����һ�ײ�ֵ�һ��
arma::mat FaceJ_dU2;//�غ�����һ�ײ�ֵڶ���
arma::mat FaceJ_dU3;//�غ�����һ�ײ�ֵ�����
arma::mat FaceJ_dU4;//�غ�����һ�ײ�ֵ�����
arma::mat FaceJ_dddU1;//�غ��������ײ�ֵ�һ��
arma::mat FaceJ_dddU2;//�غ��������ײ�ֵڶ���
arma::mat FaceJ_dddU3;//�غ��������ײ�ֵ�����
arma::mat FaceJ_dddU4;//�غ��������ײ�ֵ�����
arma::mat FaceJ_D1;//���ߵ��˹�ճ��ͨ����һ��
arma::mat FaceJ_D2;//���ߵ��˹�ճ��ͨ���ڶ���
arma::mat FaceJ_D3;//���ߵ��˹�ճ��ͨ��������
arma::mat FaceJ_D4;//���ߵ��˹�ճ��ͨ��������
arma::mat FaceJ_nx;//�ߵķ��������x����
arma::mat FaceJ_ny;//�ߵķ��������x����

//3.�����
arma::mat Node_Den;//�ܶ�
arma::mat Node_u;//x���ٶ�
arma::mat Node_v;//y���ٶ�
arma::mat Node_E;//������
arma::mat Node_P;//ѹ��
arma::mat Node_T;//�¶�
arma::mat Node_Ma;//�����
arma::mat Node_x;
arma::mat Node_y;

//4.��Ԫ
//4.1I����׶ˣ�south
arma::mat DummyIS_Den;//�ܶ�
arma::mat DummyIS_u;//x���ٶ�
arma::mat DummyIS_v;//y���ٶ�
arma::mat DummyIS_E;//������
arma::mat DummyIS_P;//ѹ��
arma::mat DummyIS_T;//�¶�
arma::mat DummyIS_Ma;//�����
arma::mat DummyIS_dt;//ʱ�䲽��
arma::mat DummyIS_c;//��������
arma::mat DummyIS_U1;//�غ�����һ��
arma::mat DummyIS_U2;//�غ����ڶ���
arma::mat DummyIS_U3;//�غ���������
arma::mat DummyIS_U4;//�غ���������
//4.2I���򶥶ˣ�north
arma::mat DummyIN_Den;//�ܶ�
arma::mat DummyIN_u;//x���ٶ�
arma::mat DummyIN_v;//y���ٶ�
arma::mat DummyIN_E;//������
arma::mat DummyIN_P;//ѹ��
arma::mat DummyIN_T;//�¶�
arma::mat DummyIN_Ma;//�����
arma::mat DummyIN_dt;//ʱ�䲽��
arma::mat DummyIN_c;//��������
arma::mat DummyIN_U1;//�غ�����һ��
arma::mat DummyIN_U2;//�غ����ڶ���
arma::mat DummyIN_U3;//�غ���������
arma::mat DummyIN_U4;//�غ���������
//4.3J������ˣ�west
arma::mat DummyJW_Den;//�ܶ�
arma::mat DummyJW_u;//x���ٶ�
arma::mat DummyJW_v;//y���ٶ�
arma::mat DummyJW_E;//������
arma::mat DummyJW_P;//ѹ��
arma::mat DummyJW_T;//�¶�
arma::mat DummyJW_Ma;//�����
arma::mat DummyJW_dt;//ʱ�䲽��
arma::mat DummyJW_c;//��������
arma::mat DummyJW_U1;//�غ�����һ��
arma::mat DummyJW_U2;//�غ����ڶ���
arma::mat DummyJW_U3;//�غ���������
arma::mat DummyJW_U4;//�غ���������
//4.4J�����Ҷˣ�east
arma::mat DummyJE_Den;//�ܶ�
arma::mat DummyJE_u;//x���ٶ�
arma::mat DummyJE_v;//y���ٶ�
arma::mat DummyJE_E;//������
arma::mat DummyJE_P;//ѹ��
arma::mat DummyJE_T;//�¶�
arma::mat DummyJE_Ma;//�����
arma::mat DummyJE_dt;//ʱ�䲽��
arma::mat DummyJE_c;//��������
arma::mat DummyJE_U1;//�غ�����һ��
arma::mat DummyJE_U2;//�غ����ڶ���
arma::mat DummyJE_U3;//�غ���������
arma::mat DummyJE_U4;//�غ���������

//������Ϣ����
int Imax, Jmax;//���������i������j�Ƿ���

//������Ϣ����
int n_dummy;//��Ԫ��ghost cell��������

//�����������Ҫ�ı���
int i, j;//�����ѭ������
double CFL;//CFL��
int n;//��ǰ������
int Iternum_max;//����������
double* RKalpha;//Runge_Kutta��ϵ��
double Resid;//�в�
arma::mat P_former;//Den_former[]���ڱ�����һ���ĸ����ܶ�,���ڲв����
arma::mat muI;//I����ļ����������ӣ��˹�ճ����
arma::mat muJ;//J����ļ����������ӣ��˹�ճ����
int numAOA;//���Ǽ���
int M;//�����

//��������ĳ���
double k;//���ȱ�
double R;//���峣��
double T0;//��ƽ���¶Ȳο�ֵ
double P0;//��ƽ��ѹ���ο�ֵ
double c0;//��ƽ�����ٲο�ֵ

//��������
double Tt;//�����¶�
double Pt;//����ѹ��
double Dent;//�����ܶ�
double ut;//����ˮƽ�ٶ�
double vt;//������ֱ�ٶ�
double Ma;//���������
double ct;//��������
double Et;//��������
double AOA;//�������ǣ���λ��
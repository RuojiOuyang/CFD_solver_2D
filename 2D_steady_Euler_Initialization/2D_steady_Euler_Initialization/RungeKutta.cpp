/*
RungeKutta������ʱ���ƽ����̣��ռ���ɢҲ����������
*/

#include "Main.h"
#include "LocalTimeStep.h"
#include "BoundaryCondition.h"
#include "Convective.h"
#include "GhostCells.h"
#include "Dissipative.h"
#include "Primitive.h"

using namespace std;

void RungeKutta()
{
	int z;//������ʹ�õ�����������5��RK�������z��1��5��

	//������һ��������
	P_former = Cell_P;
	Cell_U1_former = Cell_U1;
	Cell_U2_former = Cell_U2;
	Cell_U3_former = Cell_U3;
	Cell_U4_former = Cell_U4;

	//RK�ļ��㣬Ӧ�ü�������ʱ���
	for (z = 1; z <= 5; z++)
	{
		LocalTimeStep();//���㵱��ʱ�䲽
		BoundaryCondition();//�߽�����
		Convective();//����������嵥Ԫ��ͨ��
		GhostCells();//������Ԫ
		Dissipative();//�˹����

		Cell_Med1 = (Cell_QC1 - Cell_D1) / Cell_Vol;
		Cell_U1 = Cell_U1_former - RKalpha[z] * Cell_dt % Cell_Med1;

		Cell_Med2 = (Cell_QC2 - Cell_D2) / Cell_Vol;
		Cell_U2 = Cell_U2_former - RKalpha[z] * Cell_dt % Cell_Med2;

		Cell_Med3 = (Cell_QC3 - Cell_D3) / Cell_Vol;
		Cell_U3 = Cell_U3_former - RKalpha[z] * Cell_dt % Cell_Med3;

		Cell_Med4 = (Cell_QC4 - Cell_D4) / Cell_Vol;
		Cell_U4 = Cell_U4_former - RKalpha[z] * Cell_dt % Cell_Med4;

		Primitive();//���غ�������Ϊԭʼ������������һ������
	}
}
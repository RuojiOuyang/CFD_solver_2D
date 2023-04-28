/*
Runge_Kutta������ʱ����ɢ���̣��ռ���ɢҲ����������
*/

#include "Main.h"
#include "Min_time_step.h"
#include "Boundary_condition.h"
#include "Convective.h"
#include "Ghost_cells.h"
#include "Diffusive.h"
#include "Dissipative.h"
#include "Primitive.h"

using namespace std;

void Runge_Kutta()
{
	int z;//������ʹ�õ�����������5��RK�������z��1:5��

	//������һ��������
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			//Den_former[]���ڱ�����һ���ĸ����ܶ�,���ڲв����
			Den_former[i][j] = Cell[i][j].Den;
			//Cell[i][j].U_former���ڱ�����һ�����غ���������Runge-Kutta�ļ���
			Cell[i][j].U_former[1] = Cell[i][j].U[1];
			Cell[i][j].U_former[2] = Cell[i][j].U[2];
			Cell[i][j].U_former[3] = Cell[i][j].U[3];
			Cell[i][j].U_former[4] = Cell[i][j].U[4];
		}
	}

	//RK���ļ���
	for (z = 1; z <= 5; z++)
	{
		Min_time_step();//������Сʱ�䲽
		Boundary_condition();//�߽�����
		Ghost_cells();//������Ԫ
		Convective();//����������嵥Ԫ�Ķ�����
		Diffusive();//����������嵥Ԫ����ɢ��
		Dissipative();//�˹����

		for (i = 1; i <= Imax - 1; i++)
		{
			for (j = 1; j <= Jmax - 1; j++)
			{
				//Cell[i][j].QV[1] = Cell[i][j].QV[2] = Cell[i][j].QV[3] = Cell[i][j].QV[4] = 0;
				//Cell[i][j].D[1] = Cell[i][j].D[2] = Cell[i][j].D[3] = Cell[i][j].D[4] = 0;
				
				Cell[i][j].Med[1] = (Cell[i][j].QC[1] - Cell[i][j].QV[1] - Cell[i][j].D[1]) / Cell[i][j].Vol;
				Cell[i][j].U[1] = Cell[i][j].U_former[1] - RKalpha[z] * Cell[i][j].dt * Cell[i][j].Med[1];

				Cell[i][j].Med[2] = (Cell[i][j].QC[2] - Cell[i][j].QV[2] - Cell[i][j].D[2]) / Cell[i][j].Vol;
				Cell[i][j].U[2] = Cell[i][j].U_former[2] - RKalpha[z] * Cell[i][j].dt * Cell[i][j].Med[2];

				Cell[i][j].Med[3] = (Cell[i][j].QC[3] - Cell[i][j].QV[3] - Cell[i][j].D[3]) / Cell[i][j].Vol;
				Cell[i][j].U[3] = Cell[i][j].U_former[3] - RKalpha[z] * Cell[i][j].dt * Cell[i][j].Med[3];

				Cell[i][j].Med[4] = (Cell[i][j].QC[4] - Cell[i][j].QV[4] - Cell[i][j].D[4]) / Cell[i][j].Vol;
				Cell[i][j].U[4] = Cell[i][j].U_former[4] - RKalpha[z] * Cell[i][j].dt * Cell[i][j].Med[4];
			}
		}
		Primitive();//���غ�������Ϊԭʼ������������һ������
	}
}
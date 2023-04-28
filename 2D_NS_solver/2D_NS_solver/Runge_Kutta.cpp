/*
Runge_Kutta函数，时间离散过程，空间离散也包含在里面
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
	int z;//迭代中使用的量，这里是5步RK法，因此z是1:5的

	//保存上一步的数据
	for (i = 1; i <= Imax - 1; i++)
	{
		for (j = 1; j <= Jmax - 1; j++)
		{
			//Den_former[]用于保存上一步的各点密度,用于残差计算
			Den_former[i][j] = Cell[i][j].Den;
			//Cell[i][j].U_former用于保存上一步的守恒量，用于Runge-Kutta的计算
			Cell[i][j].U_former[1] = Cell[i][j].U[1];
			Cell[i][j].U_former[2] = Cell[i][j].U[2];
			Cell[i][j].U_former[3] = Cell[i][j].U[3];
			Cell[i][j].U_former[4] = Cell[i][j].U[4];
		}
	}

	//RK法的计算
	for (z = 1; z <= 5; z++)
	{
		Min_time_step();//计算最小时间步
		Boundary_condition();//边界条件
		Ghost_cells();//计算虚元
		Convective();//计算出控制体单元的对流项
		Diffusive();//计算出控制体单元的扩散项
		Dissipative();//人工黏性

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
		Primitive();//将守恒量解耦为原始变量，用于下一步计算
	}
}
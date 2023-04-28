/*
RungeKutta函数，时间推进过程，空间离散也包含在里面
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
	int z;//迭代中使用的量，这里是5步RK法，因此z是1：5的

	//保存上一步的数据
	P_former = Cell_P;
	Cell_U1_former = Cell_U1;
	Cell_U2_former = Cell_U2;
	Cell_U3_former = Cell_U3;
	Cell_U4_former = Cell_U4;

	//RK的计算，应该计算所有时间层
	for (z = 1; z <= 5; z++)
	{
		LocalTimeStep();//计算当地时间步
		BoundaryCondition();//边界条件
		Convective();//计算出控制体单元的通量
		GhostCells();//计算虚元
		Dissipative();//人工黏性

		Cell_Med1 = (Cell_QC1 - Cell_D1) / Cell_Vol;
		Cell_U1 = Cell_U1_former - RKalpha[z] * Cell_dt % Cell_Med1;

		Cell_Med2 = (Cell_QC2 - Cell_D2) / Cell_Vol;
		Cell_U2 = Cell_U2_former - RKalpha[z] * Cell_dt % Cell_Med2;

		Cell_Med3 = (Cell_QC3 - Cell_D3) / Cell_Vol;
		Cell_U3 = Cell_U3_former - RKalpha[z] * Cell_dt % Cell_Med3;

		Cell_Med4 = (Cell_QC4 - Cell_D4) / Cell_Vol;
		Cell_U4 = Cell_U4_former - RKalpha[z] * Cell_dt % Cell_Med4;

		Primitive();//将守恒量解耦为原始变量，用于下一步计算
	}
}
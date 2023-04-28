/*
Geometry函数用来根据读取的网格，计算控制体体积，形心，面积矢量（式(2-1-13)出现的）
控制体体积存在中心点处，面积矢量存在边界处
*/

#include "Main.h"

using namespace std;

void Geometry()
{
	double A, B;//计算过程中的中间量，局部变量
	//double a;
	//double b;
	//double x1, x2, x3, x4, y1, y2, y3, y4;

	//求控制体体积
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			A = (Node[i + 1][j + 1].x - Node[i][j].x) * (Node[i][j + 1].y - Node[i + 1][j].y);
			B = (Node[i][j + 1].x - Node[i + 1][j].x) * (Node[i + 1][j + 1].y - Node[i][j].y);
			Cell[i][j].Vol = 0.5 * (A - B);
		}
	}
	
	//求形心位置
	/*for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			x1 = Node[i][j].x;
			x2 = Node[i + 1][j].x;
			x3 = Node[i + 1][j + 1].x;
			x4 = Node[i][j + 1].x;
			y1 = Node[i][j].y;
			y2 = Node[i + 1][j].y;
			y3 = Node[i + 1][j + 1].y;
			y4 = Node[i][j + 1].y;

			a = (x1 + x2) * (x1 * y2 - x2 * y1) + (x2 + x3) * (x2 * y3 - x3 * y2) + (x3 + x4) * (x3 * y4 - x4 * y3) + (x4 + x1) * (x4 * y1 - x1 * y4);
			b = (y1 + y2) * (x1 * y2 - x2 * y1) + (y2 + y3) * (x2 * y3 - x3 * y2) + (y3 + y4) * (x3 * y4 - x4 * y3) + (y4 + y1) * (x4 * y1 - x1 * y4);

			Cell[i][j].x = a / 6.0 / Cell[i][j].Vol;
			Cell[i][j].y = b / 6.0 / Cell[i][j].Vol;
		}
	}*/

	//求面积矢量
	for (j = 1; j <= Jmax; j++)//j从1开始到Jmax，对应着中心点上下边的边界
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			FaceI[i][j].nx = Node[i][j].y - Node[i + 1][j].y;
			FaceI[i][j].ny = Node[i + 1][j].x - Node[i][j].x;
		}
	}
	/*
	以4*4的网格为例，相当于
	- - -
	- - -
	- - -
	- - -
	*/
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax; i++)//i从1开始到Imax，对应着中心点左右边的边界
		{
			FaceJ[i][j].nx = Node[i][j + 1].y - Node[i][j].y;
			FaceJ[i][j].ny = Node[i][j].x - Node[i][j + 1].x;
		}
	}
	/*
	变成了
	 _ _ _
	|_|_|_|
	|_|_|_|
	|_|_|_|
	就画出了这个网格
	*/
}
/*
Geometry函数用来根据读取的网格，计算控制体体积，形心，面积矢量
控制体体积存在中心点处，面积矢量存在边界处
*/

#include "Main.h"

using namespace std;
using namespace arma;

void Geometry()
{
	mat a, b;//计算过程中的中间量，局部变量
	mat x1, x2, x3, x4, y1, y2, y3, y4;
	mat AA, BB;

	//1.求控制体体积
	AA = (Node_x(span(1, (long long)Imax - 1), span(1, (long long)Jmax - 1)) - Node_x(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 2))) %
		(Node_y(span(0, (long long)Imax - 2), span(1, (long long)Jmax - 1)) - Node_y(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 2)));
	BB = (Node_x(span(0, (long long)Imax - 2), span(1, (long long)Jmax - 1)) - Node_x(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 2))) %
		(Node_y(span(1, (long long)Imax - 1), span(1, (long long)Jmax - 1)) - Node_y(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 2)));
	Cell_Vol = 0.5 * (AA - BB);


	//2.求形心位置
	x1 = Node_x(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 2));
	x2 = Node_x(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 2));
	x3 = Node_x(span(1, (long long)Imax - 1), span(1, (long long)Jmax - 1));
	x4 = Node_x(span(0, (long long)Imax - 2), span(1, (long long)Jmax - 1));
	y1 = Node_y(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 2));
	y2 = Node_y(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 2));
	y3 = Node_y(span(1, (long long)Imax - 1), span(1, (long long)Jmax - 1));
	y4 = Node_y(span(0, (long long)Imax - 2), span(1, (long long)Jmax - 1));

	a = (x1 + x2) % (x1 % y2 - x2 % y1) + (x2 + x3) % (x2 % y3 - x3 % y2) + (x3 + x4) % (x3 % y4 - x4 % y3) + (x4 + x1) % (x4 % y1 - x1 % y4);
	b = (y1 + y2) % (x1 % y2 - x2 % y1) + (y2 + y3) % (x2 % y3 - x3 % y2) + (y3 + y4) % (x3 % y4 - x4 % y3) + (y4 + y1) % (x4 % y1 - x1 % y4);
	Cell_x = a / 6.0 / Cell_Vol;
	Cell_y = b / 6.0 / Cell_Vol;

	//3.求面积矢量
	FaceI_nx = Node_y(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 1)) - Node_y(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 1));
	FaceI_ny = Node_x(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 1)) - Node_x(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 1));
	FaceJ_nx = Node_y(span(0, (long long)Imax - 1), span(1, (long long)Jmax - 1)) - Node_y(span(0, (long long)Imax - 1), span(0, (long long)Jmax - 2));
	FaceJ_ny = Node_x(span(0, (long long)Imax - 1), span(0, (long long)Jmax - 2)) - Node_x(span(0, (long long)Imax - 1), span(1, (long long)Jmax - 1));
}
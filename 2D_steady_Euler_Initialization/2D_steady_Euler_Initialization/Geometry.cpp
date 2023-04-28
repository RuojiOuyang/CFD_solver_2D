/*
Geometry�����������ݶ�ȡ�����񣬼����������������ģ����ʸ��
����������������ĵ㴦�����ʸ�����ڱ߽紦
*/

#include "Main.h"

using namespace std;
using namespace arma;

void Geometry()
{
	mat a, b;//��������е��м������ֲ�����
	mat x1, x2, x3, x4, y1, y2, y3, y4;
	mat AA, BB;

	//1.����������
	AA = (Node_x(span(1, (long long)Imax - 1), span(1, (long long)Jmax - 1)) - Node_x(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 2))) %
		(Node_y(span(0, (long long)Imax - 2), span(1, (long long)Jmax - 1)) - Node_y(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 2)));
	BB = (Node_x(span(0, (long long)Imax - 2), span(1, (long long)Jmax - 1)) - Node_x(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 2))) %
		(Node_y(span(1, (long long)Imax - 1), span(1, (long long)Jmax - 1)) - Node_y(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 2)));
	Cell_Vol = 0.5 * (AA - BB);


	//2.������λ��
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

	//3.�����ʸ��
	FaceI_nx = Node_y(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 1)) - Node_y(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 1));
	FaceI_ny = Node_x(span(1, (long long)Imax - 1), span(0, (long long)Jmax - 1)) - Node_x(span(0, (long long)Imax - 2), span(0, (long long)Jmax - 1));
	FaceJ_nx = Node_y(span(0, (long long)Imax - 1), span(1, (long long)Jmax - 1)) - Node_y(span(0, (long long)Imax - 1), span(0, (long long)Jmax - 2));
	FaceJ_ny = Node_x(span(0, (long long)Imax - 1), span(0, (long long)Jmax - 2)) - Node_x(span(0, (long long)Imax - 1), span(1, (long long)Jmax - 1));
}
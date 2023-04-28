/*
Geometry�����������ݶ�ȡ�����񣬼����������������ģ����ʸ����ʽ(2-1-13)���ֵģ�
����������������ĵ㴦�����ʸ�����ڱ߽紦
*/

#include "Main.h"

using namespace std;

void Geometry()
{
	double A, B;//��������е��м������ֲ�����
	//double a;
	//double b;
	//double x1, x2, x3, x4, y1, y2, y3, y4;

	//����������
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			A = (Node[i + 1][j + 1].x - Node[i][j].x) * (Node[i][j + 1].y - Node[i + 1][j].y);
			B = (Node[i][j + 1].x - Node[i + 1][j].x) * (Node[i + 1][j + 1].y - Node[i][j].y);
			Cell[i][j].Vol = 0.5 * (A - B);
		}
	}
	
	//������λ��
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

	//�����ʸ��
	for (j = 1; j <= Jmax; j++)//j��1��ʼ��Jmax����Ӧ�����ĵ����±ߵı߽�
	{
		for (i = 1; i <= Imax - 1; i++)
		{
			FaceI[i][j].nx = Node[i][j].y - Node[i + 1][j].y;
			FaceI[i][j].ny = Node[i + 1][j].x - Node[i][j].x;
		}
	}
	/*
	��4*4������Ϊ�����൱��
	- - -
	- - -
	- - -
	- - -
	*/
	for (j = 1; j <= Jmax - 1; j++)
	{
		for (i = 1; i <= Imax; i++)//i��1��ʼ��Imax����Ӧ�����ĵ����ұߵı߽�
		{
			FaceJ[i][j].nx = Node[i][j + 1].y - Node[i][j].y;
			FaceJ[i][j].ny = Node[i][j].x - Node[i][j + 1].x;
		}
	}
	/*
	�����
	 _ _ _
	|_|_|_|
	|_|_|_|
	|_|_|_|
	�ͻ������������
	*/
}
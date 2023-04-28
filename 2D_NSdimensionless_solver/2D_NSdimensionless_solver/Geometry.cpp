/*
Geometry�����������ݶ�ȡ�����񣬼����������������ģ����ʸ����ʽ(2-1-13)���ֵģ�
����������������ĵ㴦�����ʸ�����ڱ߽紦
*/

#include "Main.h"

using namespace std;

void Geometry()
{
	double A, B;//��������е��м������ֲ�����

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
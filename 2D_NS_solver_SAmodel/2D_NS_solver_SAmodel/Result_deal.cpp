/*
Result_deal函数处理数据
因为只得到了中心点处的数据，还要把它转化到节点上
*/

#include "Main.h"

using namespace std;

void Result_deal()
{
	//4个边界顶点单独处理，这是因为顶点处没设虚元
	Node[1][1].Den = 1.0 / 3.0 * (DummyIS[1][1].Den + Cell[1][1].Den + DummyJW[1][1].Den);
	Node[1][1].u = 1.0 / 3.0 * (DummyIS[1][1].u + Cell[1][1].u + DummyJW[1][1].u);
	Node[1][1].v = 1.0 / 3.0 * (DummyIS[1][1].v + Cell[1][1].v + DummyJW[1][1].v);
	Node[1][1].P = 1.0 / 3.0 * (DummyIS[1][1].P + Cell[1][1].P + DummyJW[1][1].P);
	Node[1][1].T = 1.0 / 3.0 * (DummyIS[1][1].T + Cell[1][1].T + DummyJW[1][1].T);
	Node[1][1].E = 1.0 / 3.0 * (DummyIS[1][1].E + Cell[1][1].E + DummyJW[1][1].E);
	Node[1][1].H = 1.0 / 3.0 * (DummyIS[1][1].H + Cell[1][1].H + DummyJW[1][1].H);
	Node[1][1].Ma = 1.0 / 3.0 * (DummyIS[1][1].Ma + Cell[1][1].Ma + DummyJW[1][1].Ma);
	Node[1][1].u_y = 1.0 / 3.0 * (DummyIS[1][1].u_y + Cell[1][1].u_y + DummyJW[1][1].u_y);
	Node[1][1].v_x = 1.0 / 3.0 * (DummyIS[1][1].v_x + Cell[1][1].v_x + DummyJW[1][1].v_x);

	Node[1][Jmax].Den = 1.0 / 3.0 * (DummyIN[1][1].Den + Cell[1][Jmax - 1].Den + DummyJW[1][Jmax - 1].Den);
	Node[1][Jmax].u = 1.0 / 3.0 * (DummyIN[1][1].u + Cell[1][Jmax - 1].u + DummyJW[1][Jmax - 1].u);
	Node[1][Jmax].v = 1.0 / 3.0 * (DummyIN[1][1].v + Cell[1][Jmax - 1].v + DummyJW[1][Jmax - 1].v);
	Node[1][Jmax].P = 1.0 / 3.0 * (DummyIN[1][1].P + Cell[1][Jmax - 1].P + DummyJW[1][Jmax - 1].P);
	Node[1][Jmax].T = 1.0 / 3.0 * (DummyIN[1][1].T + Cell[1][Jmax - 1].T + DummyJW[1][Jmax - 1].T);
	Node[1][Jmax].E = 1.0 / 3.0 * (DummyIN[1][1].E + Cell[1][Jmax - 1].E + DummyJW[1][Jmax - 1].E);
	Node[1][Jmax].H = 1.0 / 3.0 * (DummyIN[1][1].H + Cell[1][Jmax - 1].H + DummyJW[1][Jmax - 1].H);
	Node[1][Jmax].Ma = 1.0 / 3.0 * (DummyIN[1][1].Ma + Cell[1][Jmax - 1].Ma + DummyJW[1][Jmax - 1].Ma);
	Node[1][Jmax].u_y = 1.0 / 3.0 * (DummyIN[1][1].u_y + Cell[1][Jmax - 1].u_y + DummyJW[1][Jmax - 1].u_y);
	Node[1][Jmax].v_x = 1.0 / 3.0 * (DummyIN[1][1].v_x + Cell[1][Jmax - 1].v_x + DummyJW[1][Jmax - 1].v_x);

	Node[Imax][1].Den = 1.0 / 3.0 * (DummyIS[Imax - 1][1].Den + Cell[Imax - 1][1].Den + DummyJE[1][1].Den);
	Node[Imax][1].u = 1.0 / 3.0 * (DummyIS[Imax - 1][1].u + Cell[Imax - 1][1].u + DummyJE[1][1].u);
	Node[Imax][1].v = 1.0 / 3.0 * (DummyIS[Imax - 1][1].v + Cell[Imax - 1][1].v + DummyJE[1][1].v);
	Node[Imax][1].P = 1.0 / 3.0 * (DummyIS[Imax - 1][1].P + Cell[Imax - 1][1].P + DummyJE[1][1].P);
	Node[Imax][1].T = 1.0 / 3.0 * (DummyIS[Imax - 1][1].T + Cell[Imax - 1][1].T + DummyJE[1][1].T);
	Node[Imax][1].E = 1.0 / 3.0 * (DummyIS[Imax - 1][1].E + Cell[Imax - 1][1].E + DummyJE[1][1].E);
	Node[Imax][1].H = 1.0 / 3.0 * (DummyIS[Imax - 1][1].H + Cell[Imax - 1][1].H + DummyJE[1][1].H);
	Node[Imax][1].Ma = 1.0 / 3.0 * (DummyIS[Imax - 1][1].Ma + Cell[Imax - 1][1].Ma + DummyJE[1][1].Ma);
	Node[Imax][1].u_y = 1.0 / 3.0 * (DummyIS[Imax - 1][1].u_y + Cell[Imax - 1][1].u_y + DummyJE[1][1].u_y);
	Node[Imax][1].v_x = 1.0 / 3.0 * (DummyIS[Imax - 1][1].v_x + Cell[Imax - 1][1].v_x + DummyJE[1][1].v_x);

	Node[Imax][Jmax].Den = 1.0 / 3.0 * (DummyIN[Imax - 1][1].Den + Cell[Imax - 1][Jmax - 1].Den + DummyJE[1][Jmax - 1].Den);
	Node[Imax][Jmax].u = 1.0 / 3.0 * (DummyIN[Imax - 1][1].u + Cell[Imax - 1][Jmax - 1].u + DummyJE[1][Jmax - 1].u);
	Node[Imax][Jmax].v = 1.0 / 3.0 * (DummyIN[Imax - 1][1].v + Cell[Imax - 1][Jmax - 1].v + DummyJE[1][Jmax - 1].v);
	Node[Imax][Jmax].P = 1.0 / 3.0 * (DummyIN[Imax - 1][1].P + Cell[Imax - 1][Jmax - 1].P + DummyJE[1][Jmax - 1].P);
	Node[Imax][Jmax].T = 1.0 / 3.0 * (DummyIN[Imax - 1][1].T + Cell[Imax - 1][Jmax - 1].T + DummyJE[1][Jmax - 1].T);
	Node[Imax][Jmax].E = 1.0 / 3.0 * (DummyIN[Imax - 1][1].E + Cell[Imax - 1][Jmax - 1].E + DummyJE[1][Jmax - 1].E);
	Node[Imax][Jmax].H = 1.0 / 3.0 * (DummyIN[Imax - 1][1].H + Cell[Imax - 1][Jmax - 1].H + DummyJE[1][Jmax - 1].H);
	Node[Imax][Jmax].Ma = 1.0 / 3.0 * (DummyIN[Imax - 1][1].Ma + Cell[Imax - 1][Jmax - 1].Ma + DummyJE[1][Jmax - 1].Ma);
	Node[Imax][Jmax].u_y = 1.0 / 3.0 * (DummyIN[Imax - 1][1].u_y + Cell[Imax - 1][Jmax - 1].u_y + DummyJE[1][Jmax - 1].u_y);
	Node[Imax][Jmax].v_x = 1.0 / 3.0 * (DummyIN[Imax - 1][1].v_x + Cell[Imax - 1][Jmax - 1].v_x + DummyJE[1][Jmax - 1].v_x);

	//4个边界
	for (j = 2; j <= Jmax - 1; j++)
	{
		Node[1][j].Den = 1.0 / 4.0 * (DummyJW[1][j - 1].Den + DummyJW[1][j].Den + Cell[1][j - 1].Den + Cell[1][j].Den);
		Node[1][j].u = 1.0 / 4.0 * (DummyJW[1][j - 1].u + DummyJW[1][j].u + Cell[1][j - 1].u + Cell[1][j].u);
		Node[1][j].v = 1.0 / 4.0 * (DummyJW[1][j - 1].v + DummyJW[1][j].v + Cell[1][j - 1].v + Cell[1][j].v);
		Node[1][j].E = 1.0 / 4.0 * (DummyJW[1][j - 1].E + DummyJW[1][j].E + Cell[1][j - 1].E + Cell[1][j].E);
		Node[1][j].H = 1.0 / 4.0 * (DummyJW[1][j - 1].H + DummyJW[1][j].H + Cell[1][j - 1].H + Cell[1][j].H);
		Node[1][j].T = 1.0 / 4.0 * (DummyJW[1][j - 1].T + DummyJW[1][j].T + Cell[1][j - 1].T + Cell[1][j].T);
		Node[1][j].P = 1.0 / 4.0 * (DummyJW[1][j - 1].P + DummyJW[1][j].P + Cell[1][j - 1].P + Cell[1][j].P);
		Node[1][j].Ma = 1.0 / 4.0 * (DummyJW[1][j - 1].Ma + DummyJW[1][j].Ma + Cell[1][j - 1].Ma + Cell[1][j].Ma);
		Node[1][j].u_y = 1.0 / 4.0 * (DummyJW[1][j - 1].u_y + DummyJW[1][j].u_y + Cell[1][j - 1].u_y + Cell[1][j].u_y);
		Node[1][j].v_x = 1.0 / 4.0 * (DummyJW[1][j - 1].v_x + DummyJW[1][j].v_x + Cell[1][j - 1].v_x + Cell[1][j].v_x);

		Node[Imax][j].Den = 1.0 / 4.0 * (DummyJE[1][j - 1].Den + DummyJE[1][j].Den + Cell[Imax - 1][j - 1].Den + Cell[Imax - 1][j].Den);
		Node[Imax][j].u = 1.0 / 4.0 * (DummyJE[1][j - 1].u + DummyJE[1][j].u + Cell[Imax - 1][j - 1].u + Cell[Imax - 1][j].u);
		Node[Imax][j].v = 1.0 / 4.0 * (DummyJE[1][j - 1].v + DummyJE[1][j].v + Cell[Imax - 1][j - 1].v + Cell[Imax - 1][j].v);
		Node[Imax][j].E = 1.0 / 4.0 * (DummyJE[1][j - 1].E + DummyJE[1][j].E + Cell[Imax - 1][j - 1].E + Cell[Imax - 1][j].E);
		Node[Imax][j].H = 1.0 / 4.0 * (DummyJE[1][j - 1].H + DummyJE[1][j].H + Cell[Imax - 1][j - 1].H + Cell[Imax - 1][j].H);
		Node[Imax][j].T = 1.0 / 4.0 * (DummyJE[1][j - 1].T + DummyJE[1][j].T + Cell[Imax - 1][j - 1].T + Cell[Imax - 1][j].T);
		Node[Imax][j].P = 1.0 / 4.0 * (DummyJE[1][j - 1].P + DummyJE[1][j].P + Cell[Imax - 1][j - 1].P + Cell[Imax - 1][j].P);
		Node[Imax][j].Ma = 1.0 / 4.0 * (DummyJE[1][j - 1].Ma + DummyJE[1][j].Ma + Cell[Imax - 1][j - 1].Ma + Cell[Imax - 1][j].Ma);
		Node[Imax][j].u_y = 1.0 / 4.0 * (DummyJE[1][j - 1].u_y + DummyJE[1][j].u_y + Cell[Imax - 1][j - 1].u_y + Cell[Imax - 1][j].u_y);
		Node[Imax][j].v_x = 1.0 / 4.0 * (DummyJE[1][j - 1].v_x + DummyJE[1][j].v_x + Cell[Imax - 1][j - 1].v_x + Cell[Imax - 1][j].v_x);
	}
	for (i = 2; i <= Imax - 1; i++)
	{
		Node[i][1].Den = 1.0 / 4.0 * (DummyIS[i - 1][1].Den + DummyIS[i][1].Den + Cell[i - 1][1].Den + Cell[i][1].Den);
		Node[i][1].u = 1.0 / 4.0 * (DummyIS[i - 1][1].u + DummyIS[i][1].u + Cell[i - 1][1].u + Cell[i][1].u);
		Node[i][1].v = 1.0 / 4.0 * (DummyIS[i - 1][1].v + DummyIS[i][1].v + Cell[i - 1][1].v + Cell[i][1].v);
		Node[i][1].E = 1.0 / 4.0 * (DummyIS[i - 1][1].E + DummyIS[i][1].E + Cell[i - 1][1].E + Cell[i][1].E);
		Node[i][1].H = 1.0 / 4.0 * (DummyIS[i - 1][1].H + DummyIS[i][1].H + Cell[i - 1][1].H + Cell[i][1].H);
		Node[i][1].T = 1.0 / 4.0 * (DummyIS[i - 1][1].T + DummyIS[i][1].T + Cell[i - 1][1].T + Cell[i][1].T);
		Node[i][1].P = 1.0 / 4.0 * (DummyIS[i - 1][1].P + DummyIS[i][1].P + Cell[i - 1][1].P + Cell[i][1].P);
		Node[i][1].Ma = 1.0 / 4.0 * (DummyIS[i - 1][1].Ma + DummyIS[i][1].Ma + Cell[i - 1][1].Ma + Cell[i][1].Ma);
		Node[i][1].u_y = 1.0 / 4.0 * (DummyIS[i - 1][1].u_y + DummyIS[i][1].u_y + Cell[i - 1][1].u_y + Cell[i][1].u_y);
		Node[i][1].v_x = 1.0 / 4.0 * (DummyIS[i - 1][1].v_x + DummyIS[i][1].v_x + Cell[i - 1][1].v_x + Cell[i][1].v_x);

		Node[i][Jmax].Den = 1.0 / 4.0 * (DummyIN[i - 1][1].Den + DummyIN[i][1].Den + Cell[i - 1][Jmax - 1].Den + Cell[i][Jmax - 1].Den);
		Node[i][Jmax].u = 1.0 / 4.0 * (DummyIN[i - 1][1].u + DummyIN[i][1].u + Cell[i - 1][Jmax - 1].u + Cell[i][Jmax - 1].u);
		Node[i][Jmax].v = 1.0 / 4.0 * (DummyIN[i - 1][1].v + DummyIN[i][1].v + Cell[i - 1][Jmax - 1].v + Cell[i][Jmax - 1].v);
		Node[i][Jmax].E = 1.0 / 4.0 * (DummyIN[i - 1][1].E + DummyIN[i][1].E + Cell[i - 1][Jmax - 1].E + Cell[i][Jmax - 1].E);
		Node[i][Jmax].H = 1.0 / 4.0 * (DummyIN[i - 1][1].H + DummyIN[i][1].H + Cell[i - 1][Jmax - 1].H + Cell[i][Jmax - 1].H);
		Node[i][Jmax].T = 1.0 / 4.0 * (DummyIN[i - 1][1].T + DummyIN[i][1].T + Cell[i - 1][Jmax - 1].T + Cell[i][Jmax - 1].T);
		Node[i][Jmax].P = 1.0 / 4.0 * (DummyIN[i - 1][1].P + DummyIN[i][1].P + Cell[i - 1][Jmax - 1].P + Cell[i][Jmax - 1].P);
		Node[i][Jmax].Ma = 1.0 / 4.0 * (DummyIN[i - 1][1].Ma + DummyIN[i][1].Ma + Cell[i - 1][Jmax - 1].Ma + Cell[i][Jmax - 1].Ma);
		Node[i][Jmax].u_y = 1.0 / 4.0 * (DummyIN[i - 1][1].u_y + DummyIN[i][1].u_y + Cell[i - 1][Jmax - 1].u_y + Cell[i][Jmax - 1].u_y);
		Node[i][Jmax].v_x = 1.0 / 4.0 * (DummyIN[i - 1][1].v_x + DummyIN[i][1].v_x + Cell[i - 1][Jmax - 1].v_x + Cell[i][Jmax - 1].v_x);
	}

	//内部
	for (i = 2; i <= Imax - 1; i++)
	{
		for (j = 2; j <= Jmax - 1; j++)
		{
			Node[i][j].Den = 1.0 / 4.0 * (Cell[i - 1][j - 1].Den + Cell[i][j - 1].Den + Cell[i - 1][j].Den + Cell[i][j].Den);
			Node[i][j].u = 1.0 / 4.0 * (Cell[i - 1][j - 1].u + Cell[i][j - 1].u + Cell[i - 1][j].u + Cell[i][j].u);
			Node[i][j].v = 1.0 / 4.0 * (Cell[i - 1][j - 1].v + Cell[i][j - 1].v + Cell[i - 1][j].v + Cell[i][j].v);
			Node[i][j].P = 1.0 / 4.0 * (Cell[i - 1][j - 1].P + Cell[i][j - 1].P + Cell[i - 1][j].P + Cell[i][j].P);
			Node[i][j].T = 1.0 / 4.0 * (Cell[i - 1][j - 1].T + Cell[i][j - 1].T + Cell[i - 1][j].T + Cell[i][j].T);
			Node[i][j].E = 1.0 / 4.0 * (Cell[i - 1][j - 1].E + Cell[i][j - 1].E + Cell[i - 1][j].E + Cell[i][j].E);
			Node[i][j].H = 1.0 / 4.0 * (Cell[i - 1][j - 1].H + Cell[i][j - 1].H + Cell[i - 1][j].H + Cell[i][j].H);
			Node[i][j].Ma = 1.0 / 4.0 * (Cell[i - 1][j - 1].Ma + Cell[i][j - 1].Ma + Cell[i - 1][j].Ma + Cell[i][j].Ma);
			Node[i][j].u_y = 1.0 / 4.0 * (Cell[i - 1][j - 1].u_y + Cell[i][j - 1].u_y + Cell[i - 1][j].u_y + Cell[i][j].u_y);
			Node[i][j].v_x = 1.0 / 4.0 * (Cell[i - 1][j - 1].v_x + Cell[i][j - 1].v_x + Cell[i - 1][j].v_x + Cell[i][j].v_x);
		}
	}
}
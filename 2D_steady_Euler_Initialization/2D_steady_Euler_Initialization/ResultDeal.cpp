/*
ResultDeal函数处理数据
因为只得到了中心点处的数据，还要把它转化到节点上
*/

#include "Main.h"

using namespace std;
using namespace arma;

void ResultDeal()
{
	//4个边界顶点单独处理，这是因为顶点处没设虚元
	Node_Den(0, 0) = 1.0 / 3.0 * (DummyIS_Den(0, 0) + Cell_Den(0, 0) + DummyJW_Den(0, 0));
	Node_u(0, 0) = 1.0 / 3.0 * (DummyIS_u(0, 0) + Cell_u(0, 0) + DummyJW_u(0, 0));
	Node_v(0, 0) = 1.0 / 3.0 * (DummyIS_v(0, 0) + Cell_v(0, 0) + DummyJW_v(0, 0));
	Node_E(0, 0) = 1.0 / 3.0 * (DummyIS_E(0, 0) + Cell_E(0, 0) + DummyJW_E(0, 0));
	Node_P(0, 0) = 1.0 / 3.0 * (DummyIS_P(0, 0) + Cell_P(0, 0) + DummyJW_P(0, 0));
	Node_T(0, 0) = 1.0 / 3.0 * (DummyIS_T(0, 0) + Cell_T(0, 0) + DummyJW_T(0, 0));
	Node_Ma(0, 0) = 1.0 / 3.0 * (DummyIS_Ma(0, 0) + Cell_Ma(0, 0) + DummyJW_Ma(0, 0));

	Node_Den(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_Den(0, 0) + Cell_Den(0, (long long)Jmax - 2) + DummyJW_Den(0, (long long)Jmax - 2));
	Node_u(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_u(0, 0) + Cell_u(0, (long long)Jmax - 2) + DummyJW_u(0, (long long)Jmax - 2));
	Node_v(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_v(0, 0) + Cell_v(0, (long long)Jmax - 2) + DummyJW_v(0, (long long)Jmax - 2));
	Node_E(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_E(0, 0) + Cell_E(0, (long long)Jmax - 2) + DummyJW_E(0, (long long)Jmax - 2));
	Node_P(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_P(0, 0) + Cell_P(0, (long long)Jmax - 2) + DummyJW_P(0, (long long)Jmax - 2));
	Node_T(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_T(0, 0) + Cell_T(0, (long long)Jmax - 2) + DummyJW_T(0, (long long)Jmax - 2));
	Node_Ma(0, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_Ma(0, 0) + Cell_Ma(0, (long long)Jmax - 2) + DummyJW_Ma(0, (long long)Jmax - 2));

	Node_Den((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_Den((long long)Imax - 2, 0) + Cell_Den((long long)Imax - 2, 0) + DummyJE_Den(0, 0));
	Node_u((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_u((long long)Imax - 2, 0) + Cell_u((long long)Imax - 2, 0) + DummyJE_u(0, 0));
	Node_v((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_v((long long)Imax - 2, 0) + Cell_v((long long)Imax - 2, 0) + DummyJE_v(0, 0));
	Node_E((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_E((long long)Imax - 2, 0) + Cell_E((long long)Imax - 2, 0) + DummyJE_E(0, 0));
	Node_P((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_P((long long)Imax - 2, 0) + Cell_P((long long)Imax - 2, 0) + DummyJE_P(0, 0));
	Node_T((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_T((long long)Imax - 2, 0) + Cell_T((long long)Imax - 2, 0) + DummyJE_T(0, 0));
	Node_Ma((long long)Imax - 1, 0) = 1.0 / 3.0 * (DummyIS_Ma((long long)Imax - 2, 0) + Cell_Ma((long long)Imax - 2, 0) + DummyJE_Ma(0, 0));

	Node_Den((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_Den((long long)Imax - 2, 0) + Cell_Den((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_Den(0, (long long)Jmax - 2));
	Node_u((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_u((long long)Imax - 2, 0) + Cell_u((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_u(0, (long long)Jmax - 2));
	Node_v((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_v((long long)Imax - 2, 0) + Cell_v((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_v(0, (long long)Jmax - 2));
	Node_E((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_E((long long)Imax - 2, 0) + Cell_E((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_E(0, (long long)Jmax - 2));
	Node_P((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_P((long long)Imax - 2, 0) + Cell_P((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_P(0, (long long)Jmax - 2));
	Node_T((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_T((long long)Imax - 2, 0) + Cell_T((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_T(0, (long long)Jmax - 2));
	Node_Ma((long long)Imax - 1, (long long)Jmax - 1) = 1.0 / 3.0 * (DummyIN_Ma((long long)Imax - 2, 0) + Cell_Ma((long long)Imax - 2, (long long)Jmax - 2) + DummyJE_Ma(0, (long long)Jmax - 2));

	//4个边界
	Node_Den(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_Den(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_Den(span(0), span(1, (long long)Jmax - 2)) + Cell_Den(span(0), span(0, (long long)Jmax - 3)) + Cell_Den(span(0), span(1, (long long)Jmax - 2)));
	Node_u(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_u(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_u(span(0), span(1, (long long)Jmax - 2)) + Cell_u(span(0), span(0, (long long)Jmax - 3)) + Cell_u(span(0), span(1, (long long)Jmax - 2)));
	Node_v(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_v(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_v(span(0), span(1, (long long)Jmax - 2)) + Cell_v(span(0), span(0, (long long)Jmax - 3)) + Cell_v(span(0), span(1, (long long)Jmax - 2)));
	Node_E(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_E(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_E(span(0), span(1, (long long)Jmax - 2)) + Cell_E(span(0), span(0, (long long)Jmax - 3)) + Cell_E(span(0), span(1, (long long)Jmax - 2)));
	Node_P(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_P(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_P(span(0), span(1, (long long)Jmax - 2)) + Cell_P(span(0), span(0, (long long)Jmax - 3)) + Cell_P(span(0), span(1, (long long)Jmax - 2)));
	Node_T(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_T(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_T(span(0), span(1, (long long)Jmax - 2)) + Cell_T(span(0), span(0, (long long)Jmax - 3)) + Cell_T(span(0), span(1, (long long)Jmax - 2)));
	Node_Ma(span(0), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJW_Ma(span(0), span(0, (long long)Jmax - 3)) +
		DummyJW_Ma(span(0), span(1, (long long)Jmax - 2)) + Cell_Ma(span(0), span(0, (long long)Jmax - 3)) + Cell_Ma(span(0), span(1, (long long)Jmax - 2)));

	Node_Den(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_Den(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_Den(span(0), span(1, (long long)Jmax - 2)) + Cell_Den(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_Den(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_u(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_u(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_u(span(0), span(1, (long long)Jmax - 2)) + Cell_u(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_u(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_v(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_v(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_v(span(0), span(1, (long long)Jmax - 2)) + Cell_v(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_v(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_E(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_E(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_E(span(0), span(1, (long long)Jmax - 2)) + Cell_E(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_E(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_P(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_P(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_P(span(0), span(1, (long long)Jmax - 2)) + Cell_P(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_P(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_T(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_T(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_T(span(0), span(1, (long long)Jmax - 2)) + Cell_T(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_T(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_Ma(span((long long)Imax - 1), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (DummyJE_Ma(span(0), span(0, (long long)Jmax - 3)) +
		DummyJE_Ma(span(0), span(1, (long long)Jmax - 2)) + Cell_Ma(span((long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_Ma(span((long long)Imax - 2), span(1, (long long)Jmax - 2)));

	Node_Den(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_Den(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_Den(span(1, (long long)Imax - 2), span(0)) + Cell_Den(span(0, (long long)Imax - 3), span(0)) + Cell_Den(span(1, (long long)Imax - 2), span(0)));
	Node_u(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_u(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_u(span(1, (long long)Imax - 2), span(0)) + Cell_u(span(0, (long long)Imax - 3), span(0)) + Cell_u(span(1, (long long)Imax - 2), span(0)));
	Node_v(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_v(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_v(span(1, (long long)Imax - 2), span(0)) + Cell_v(span(0, (long long)Imax - 3), span(0)) + Cell_v(span(1, (long long)Imax - 2), span(0)));
	Node_E(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_E(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_E(span(1, (long long)Imax - 2), span(0)) + Cell_E(span(0, (long long)Imax - 3), span(0)) + Cell_E(span(1, (long long)Imax - 2), span(0)));
	Node_P(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_P(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_P(span(1, (long long)Imax - 2), span(0)) + Cell_P(span(0, (long long)Imax - 3), span(0)) + Cell_P(span(1, (long long)Imax - 2), span(0)));
	Node_T(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_T(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_T(span(1, (long long)Imax - 2), span(0)) + Cell_T(span(0, (long long)Imax - 3), span(0)) + Cell_T(span(1, (long long)Imax - 2), span(0)));
	Node_Ma(span(1, (long long)Imax - 2), span(0)) = 1.0 / 4.0 * (DummyIS_Ma(span(0, (long long)Imax - 3), span(0)) +
		DummyIS_Ma(span(1, (long long)Imax - 2), span(0)) + Cell_Ma(span(0, (long long)Imax - 3), span(0)) + Cell_Ma(span(1, (long long)Imax - 2), span(0)));

	Node_Den(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_Den(span(0, (long long)Imax - 3), span(0)) + DummyIN_Den(span(1, (long long)Imax - 2), span(0)) +
		Cell_Den(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_Den(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));
	Node_u(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_u(span(0, (long long)Imax - 3), span(0)) + DummyIN_u(span(1, (long long)Imax - 2), span(0)) +
		Cell_u(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_u(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));
	Node_v(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_v(span(0, (long long)Imax - 3), span(0)) + DummyIN_v(span(1, (long long)Imax - 2), span(0)) +
		Cell_v(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_v(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));
	Node_E(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_E(span(0, (long long)Imax - 3), span(0)) + DummyIN_E(span(1, (long long)Imax - 2), span(0)) +
		Cell_E(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_E(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));
	Node_P(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_P(span(0, (long long)Imax - 3), span(0)) + DummyIN_P(span(1, (long long)Imax - 2), span(0)) +
		Cell_P(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_P(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));
	Node_T(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_T(span(0, (long long)Imax - 3), span(0)) + DummyIN_T(span(1, (long long)Imax - 2), span(0)) +
		Cell_T(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_T(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));
	Node_Ma(span(1, (long long)Imax - 2), span((long long)Jmax - 1)) = 1.0 / 4.0 * (DummyIN_Ma(span(0, (long long)Imax - 3), span(0)) + DummyIN_Ma(span(1, (long long)Imax - 2), span(0)) +
		Cell_Ma(span(0, (long long)Imax - 3), span((long long)Jmax - 2)) + Cell_Ma(span(1, (long long)Imax - 2), span((long long)Jmax - 2)));

	//内部
	Node_Den(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_Den(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_Den(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_Den(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_Den(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_u(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_u(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_u(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_u(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_u(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_v(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_v(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_v(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_v(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_v(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_E(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_E(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_E(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_E(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_E(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_P(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_P(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_P(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_P(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_P(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_T(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_T(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_T(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_T(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_T(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
	Node_Ma(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)) = 1.0 / 4.0 * (Cell_Ma(span(0, (long long)Imax - 3), span(0, (long long)Jmax - 3)) +
		Cell_Ma(span(1, (long long)Imax - 2), span(0, (long long)Jmax - 3)) + Cell_Ma(span(0, (long long)Imax - 3), span(1, (long long)Jmax - 2)) +
		Cell_Ma(span(1, (long long)Imax - 2), span(1, (long long)Jmax - 2)));
}
/*
Initialization������ʼ��������M��ʱ��㣩�������ñ������غ���
*/

#include "Main.h"

using namespace std;
using namespace arma;

void Initialization()
{
	//������ֵ������ʼ�����������������ֵ��Ϊ����ֵ
	Cell_Den.fill(Dent);
	Cell_u.fill(ut);
	Cell_v.fill(vt);
	Cell_E.fill(Et);
	Cell_P.fill(Pt);
	Cell_T.fill(Tt);
	Cell_c.fill(ct);
	Cell_Ma.fill(sqrt(ut * ut + vt * ct) / ct);

	Cell_U1.fill(Dent);
	Cell_U2.fill(Dent * ut);
	Cell_U3.fill(Dent * vt);
	Cell_U4.fill(Dent * Et);
}
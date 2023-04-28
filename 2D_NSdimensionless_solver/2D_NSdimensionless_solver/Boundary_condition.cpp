/*
Boundary_condition函数给出边界处的值，边界上的通量Face[][].FF，用在Convective里面计算Qij
*/

#include "Main.h"

using namespace std;

void Boundary_condition()
{
	//1.远场边界条件
	//i=1;Imax-1,j=Jmax为远场，远场是对Face而言的
	double c, qtinf, qninf, qti, qni, qt, qn;
	double norm_S;//方向矢量的模
	double qn_if, a_if;//判断进出口与亚超声速的量
	double Re_z, Rinf_f;

	for (i = 1; i <= Imax - 1; i++)
	{
		norm_S = sqrt(FaceI[i][Jmax].nx * FaceI[i][Jmax].nx + FaceI[i][Jmax].ny * FaceI[i][Jmax].ny);//方向矢量的模

		qtinf = (utinf * FaceI[i][Jmax].ny - vtinf * FaceI[i][Jmax].nx) / norm_S;//来流的切向速度
		qninf = (utinf * FaceI[i][Jmax].nx + vtinf * FaceI[i][Jmax].ny) / norm_S;//来流的法向速度
		qti = (Cell[i][Jmax - 1].u * ctinf * FaceI[i][Jmax].ny - Cell[i][Jmax - 1].v * ctinf * FaceI[i][Jmax].nx) / norm_S;//内场的切向速度
		qni = (Cell[i][Jmax - 1].u * ctinf * FaceI[i][Jmax].nx + Cell[i][Jmax - 1].v * ctinf * FaceI[i][Jmax].ny) / norm_S;//内场的法向速度

		if (sqrt(ut * ut + vt * vt) < ct)//亚声速
		{
			Re_z = qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1.0);
			Rinf_f = qninf - 2 * ct * ctinf / (k - 1.0);

			qn_if = 1.0 / 2.0 * (Re_z + Rinf_f);
			a_if = (k - 1.0) / 4 * (Re_z - Rinf_f);

			if (qn_if <= 0)//亚声速进口
			{
				//熵取来流
				//切向速度取来流
				qt = qtinf;
				//法向速度由内外场决定
				qn = 1.0 / 2.0 * (Re_z + Rinf_f);
				//声速由内外场决定
				c = (k - 1.0) / 4 * (Re_z - Rinf_f);
				//利用上面的量算出边界上的Den、u、v、E
				FaceI[i][Jmax].Den = pow(pow(Dentinf, k) * c * c / k / Ptinf, 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R / Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf/ ctinf;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
				//cout << FaceI[i][Jmax].E << endl;
			}
			else//亚声速出口
			{
				//熵取内场
				//切向速度取内场
				qt = qti;
				//法向速度由内外场决定
				qn = 1.0 / 2.0 * (Re_z + Rinf_f);
				//声速由内外场决定
				c = (k - 1.0) / 4 * (Re_z - Rinf_f);
				//利用上面的量算出边界上的Den、u、v、E
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den * Dentinf, k) * c * c / k / (Cell[i][Jmax - 1].P * Dentinf * ctinf * ctinf), 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R / Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf / ctinf;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
		}
		else//超声速
		{
			if (ut * FaceI[i][Jmax].nx + vt * FaceI[i][Jmax].ny < 0)//超声速进口
			{
				//熵取来流
				//切向速度取来流
				qt = qtinf;
				//法向速度由来流决定
				qn = 1.0 / 2.0 * (qninf + 2 * ctinf / (k - 1.0) + qninf - 2 * ctinf / (k - 1.0));
				//声速由来流决定
				c = (k - 1.0) / 4 * (qninf + 2 * ctinf / (k - 1.0) - qninf + 2 * ctinf / (k - 1.0));
				//利用上面的量算出边界上的Den、u、v、E
				FaceI[i][Jmax].Den = pow(pow(Dentinf, k) * c * c / k / Ptinf, 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R /Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf / ctinf;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
			else//超声速出口
			{
				//熵取内场
				//切向速度取内场
				qt = qti;
				//法向速度由内场决定
				qn = 1.0 / 2.0 * (qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1.0) + qni - 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1.0));
				//声速由内场决定
				c = (k - 1) / 4 * (qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1) - qni + 2 * Cell[i][Jmax - 1].c * ctinf / (k - 1));
				//利用上面的量算出边界上的Den、u、v、E
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den * Dentinf, k) * c * c / k / (Cell[i][Jmax - 1].P * Dentinf * ctinf * ctinf), 1.0 / (k - 1.0)) / Dentinf;
				FaceI[i][Jmax].T = c * c / k / R / Ttinf;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T * Ttinf / ctinf / ctinf;;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S / ctinf;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S / ctinf;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
			}
		}
	}
	/*
	//1.2远场边界计算守恒量
	for (i = 1; i <= Imax - 1; i++)
	{
		FaceI[i][Jmax].FU[1] = FaceI[i][Jmax].Den;
		FaceI[i][Jmax].FU[2] = FaceI[i][Jmax].Den * FaceI[i][Jmax].u;
		FaceI[i][Jmax].FU[3] = FaceI[i][Jmax].Den * FaceI[i][Jmax].v;
		FaceI[i][Jmax].FU[4] = FaceI[i][Jmax].Den * FaceI[i][Jmax].E;
	}
	//1.3由远场边界的守恒量计算通量
	double A, B;
	for (i = 1; i <= Imax - 1; i++)
	{
		A = FaceI[i][Jmax].FU[1] * (pow(FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1], 2.0) + pow(FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1], 2.0)) / 2.0;//A是den*动能
		B = (k - 1.0) * (FaceI[i][Jmax].FU[4] - A);//B是压力
		FaceI[i][Jmax].FF[1] = FaceI[i][Jmax].FU[1] * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny);
		FaceI[i][Jmax].FF[2] = FaceI[i][Jmax].FU[2] * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny) + B * FaceI[i][Jmax].nx;
		FaceI[i][Jmax].FF[3] = FaceI[i][Jmax].FU[3] * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny) + B * FaceI[i][Jmax].ny;
		FaceI[i][Jmax].FF[4] = (FaceI[i][Jmax].FU[4] + B) * (FaceI[i][Jmax].FU[2] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].nx + FaceI[i][Jmax].FU[3] / FaceI[i][Jmax].FU[1] * FaceI[i][Jmax].ny);
	}
	*/
}
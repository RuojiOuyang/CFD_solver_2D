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
	double S1, S2, S3;

	for (i = 1; i <= Imax - 1; i++)
	{
		norm_S = sqrt(FaceI[i][Jmax].nx * FaceI[i][Jmax].nx + FaceI[i][Jmax].ny * FaceI[i][Jmax].ny);//方向矢量的模

		qtinf = (ut * FaceI[i][Jmax].ny - vt * FaceI[i][Jmax].nx) / norm_S;//来流的切向速度
		qninf = (ut * FaceI[i][Jmax].nx + vt * FaceI[i][Jmax].ny) / norm_S;//来流的法向速度
		qti = (Cell[i][Jmax - 1].u * FaceI[i][Jmax].ny - Cell[i][Jmax - 1].v * FaceI[i][Jmax].nx) / norm_S;//内场的切向速度
		qni = (Cell[i][Jmax - 1].u * FaceI[i][Jmax].nx + Cell[i][Jmax - 1].v * FaceI[i][Jmax].ny) / norm_S;//内场的法向速度

		if (sqrt(ut * ut + vt * vt) < ct)//亚声速
		{
			Re_z = qni + 2 * Cell[i][Jmax - 1].c / (k - 1.0);
			Rinf_f = qninf - 2 * ct / (k - 1.0);

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
				FaceI[i][Jmax].Den = pow(pow(Dent, k) * c * c / k / Pt, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
				FaceI[i][Jmax].miubl = 0;//来流的miu一波浪等于0
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
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den, k) * c * c / k / Cell[i][Jmax - 1].P, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;

				S3 = sqrt(pow((Cell[i][Jmax - 3].x - 0.5 * Node[i][Jmax].x - 0.5 * Node[i + 1][Jmax].x), 2) + pow((Cell[i][Jmax - 3].y - 0.5 * Node[i][Jmax].y - 0.5 * Node[i + 1][Jmax].y), 2));
				S2 = sqrt(pow((Cell[i][Jmax - 2].x - 0.5 * Node[i][Jmax].x - 0.5 * Node[i + 1][Jmax].x), 2) + pow((Cell[i][Jmax - 2].y - 0.5 * Node[i][Jmax].y - 0.5 * Node[i + 1][Jmax].y), 2));
				S1 = sqrt(pow((Cell[i][Jmax - 1].x - 0.5 * Node[i][Jmax].x - 0.5 * Node[i + 1][Jmax].x), 2) + pow((Cell[i][Jmax - 1].y - 0.5 * Node[i][Jmax].y - 0.5 * Node[i + 1][Jmax].y), 2));
				FaceI[i][Jmax].miubl = (S3 / (S3 - S1) * S2 / (S2 - S1)) * Cell[i][Jmax - 1].miubl - (S3 / (S3 - S2) * S1 / (S2 - S1)) * Cell[i][Jmax - 2].miubl + (S2 / (S3 - S2) * S1 / (S3 - S1)) * Cell[i][Jmax - 3].miubl;
				if (FaceI[i][Jmax].miubl < 0)
				{
					FaceI[i][Jmax].miubl = 1e-10;
				}
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
				qn = 1.0 / 2.0 * (qninf + 2 * ct / (k - 1.0) + qninf - 2 * ct / (k - 1.0));
				//声速由来流决定
				c = (k - 1.0) / 4 * (qninf + 2 * ct / (k - 1.0) - qninf + 2 * ct / (k - 1.0));
				//利用上面的量算出边界上的Den、u、v、E
				FaceI[i][Jmax].Den = pow(pow(Dent, k) * c * c / k / Pt, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;
				FaceI[i][Jmax].miubl = 0;//来流的miu一波浪
			}
			else//超声速出口
			{
				//熵取内场
				//切向速度取内场
				qt = qti;
				//法向速度由内场决定
				qn = 1.0 / 2.0 * (qni + 2 * Cell[i][Jmax - 1].c / (k - 1.0) + qni - 2 * Cell[i][Jmax - 1].c / (k - 1.0));
				//声速由内场决定
				c = (k - 1) / 4 * (qni + 2 * Cell[i][Jmax - 1].c / (k - 1) - qni + 2 * Cell[i][Jmax - 1].c / (k - 1));
				//利用上面的量算出边界上的Den、u、v、E
				FaceI[i][Jmax].Den = pow(pow(Cell[i][Jmax - 1].Den, k) * c * c / k / Cell[i][Jmax - 1].P, 1.0 / (k - 1.0));
				FaceI[i][Jmax].T = c * c / k / R;
				FaceI[i][Jmax].P = FaceI[i][Jmax].Den * R * FaceI[i][Jmax].T;
				FaceI[i][Jmax].u = (qt * FaceI[i][Jmax].ny + qn * FaceI[i][Jmax].nx) / norm_S;
				FaceI[i][Jmax].v = (-qt * FaceI[i][Jmax].nx + qn * FaceI[i][Jmax].ny) / norm_S;
				FaceI[i][Jmax].E = FaceI[i][Jmax].P / (k - 1.0) / FaceI[i][Jmax].Den + (FaceI[i][Jmax].u * FaceI[i][Jmax].u + FaceI[i][Jmax].v * FaceI[i][Jmax].v) / 2.0;

				S3 = sqrt(pow((Cell[i][Jmax - 3].x - 0.5 * Node[i][Jmax].x - 0.5 * Node[i + 1][Jmax].x), 2) + pow((Cell[i][Jmax - 3].y - 0.5 * Node[i][Jmax].y - 0.5 * Node[i + 1][Jmax].y), 2));
				S2 = sqrt(pow((Cell[i][Jmax - 2].x - 0.5 * Node[i][Jmax].x - 0.5 * Node[i + 1][Jmax].x), 2) + pow((Cell[i][Jmax - 2].y - 0.5 * Node[i][Jmax].y - 0.5 * Node[i + 1][Jmax].y), 2));
				S1 = sqrt(pow((Cell[i][Jmax - 1].x - 0.5 * Node[i][Jmax].x - 0.5 * Node[i + 1][Jmax].x), 2) + pow((Cell[i][Jmax - 1].y - 0.5 * Node[i][Jmax].y - 0.5 * Node[i + 1][Jmax].y), 2));
				FaceI[i][Jmax].miubl = (S3 / (S3 - S1) * S2 / (S2 - S1)) * Cell[i][Jmax - 1].miubl - (S3 / (S3 - S2) * S1 / (S2 - S1)) * Cell[i][Jmax - 2].miubl + (S2 / (S3 - S2) * S1 / (S3 - S1)) * Cell[i][Jmax - 3].miubl;
				if (FaceI[i][Jmax].miubl < 0)
				{
					FaceI[i][Jmax].miubl = 1e-10;
				}
			}
		}
	}
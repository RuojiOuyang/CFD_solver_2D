/*
LocalTimeStep�������㵱��ʱ�䲽������Ϊʱ���
*/

#include "Main.h"

using namespace std;
using namespace arma;

void LocalTimeStep()
{
	//���㵱��ʱ�䲽ʱ���м����
	mat AA, BB, CC, DD;
	mat EE, FF, GG, LL, OO, PP, RR;

	AA = 0.5 * (FaceJ_nx.rows(0, (long long)Imax - 2) + FaceJ_nx.rows(1, (long long)Imax - 1));//���ʽ�Ӿ��ǰ�ʽ(2-1-13)�����һ���SI���x����ļ���һ����
	BB = 0.5 * (FaceJ_ny.rows(0, (long long)Imax - 2) + FaceJ_ny.rows(1, (long long)Imax - 1));//���ʽ�Ӿ��ǰ�ʽ(2-1-13)�����һ���SI���y����ļ���һ����
	CC = 0.5 * (FaceI_nx.cols(0, (long long)Jmax - 2) + FaceI_nx.cols(1, (long long)Jmax - 1));//���ʽ�Ӿ��ǰ�ʽ(2-1-13)����ڶ����SI���x����ļ���һ����
	DD = 0.5 * (FaceI_ny.cols(0, (long long)Jmax - 2) + FaceI_ny.cols(1, (long long)Jmax - 1));//���ʽ�Ӿ��ǰ�ʽ(2-1-13)����ڶ����SI���y����ļ���һ����

	EE = arma::abs(Cell_u % AA + Cell_v % BB);//ʽ(2-1-13)�����һ��
	FF = arma::abs(Cell_u % CC + Cell_v % DD);//ʽ(2-1-13)����ڶ���
	GG = arma::abs(arma::sqrt(AA % AA + BB % BB));
	LL = arma::abs(arma::sqrt(CC % CC + DD % DD));
	OO = EE + Cell_c % GG;
	PP = FF + Cell_c % LL;

	RR = OO + PP;

	Cell_dt = CFL * Cell_Vol / RR;//���������ʱ�䲽���������ĵ㴦
}
#include "Main.h"
#include "InitialValues.h"
#include "Initialization.h"
#include "RungeKutta.h"
#include "Residual.h"
#include "ResultDeal.h"
#include "PutOut.h"
#include "DeleteArray.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    InitialValues();//������ʼ״̬
    Initialization();//��ʼ�������������ñ������غ���

    for (n = 1; n <= Iternum_max; n++)
    {
        RungeKutta();//RK����ǰ���
        Residual();//����в�
        cout << n << "   " << Resid << endl;
        if (Resid <= 1e-3)
        {
            break;
        }
    }

    ResultDeal();//��������

    PutOut();//�������

    for (n = Iternum_max + 1;; n++)
    {
        RungeKutta();//RK����ǰ���
        Residual();//����в�
        cout << n << "   " << Resid << endl;
        if (Resid <= 1e-3)
        {
            ResultDeal();//��������
            PutOut();//�������
            break;
        }
        /*if (n % 10000 == 0)
        {
            ResultDeal();//��������
            PutOut();//�������
        }
        if (n >= 40000)
        {
            break;
        }*/
    }

    DeleteArray();//�ͷ��ڴ�
    return 0;
}
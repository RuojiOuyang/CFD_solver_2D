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
    InitialValues();//给出初始状态
    Initialization();//初始化流场，并利用变量求守恒量

    for (n = 1; n <= Iternum_max; n++)
    {
        RungeKutta();//RK法向前求解
        Residual();//计算残差
        cout << n << "   " << Resid << endl;
        if (Resid <= 1e-3)
        {
            break;
        }
    }

    ResultDeal();//处理数据

    PutOut();//输出数据

    for (n = Iternum_max + 1;; n++)
    {
        RungeKutta();//RK法向前求解
        Residual();//计算残差
        cout << n << "   " << Resid << endl;
        if (Resid <= 1e-3)
        {
            ResultDeal();//处理数据
            PutOut();//输出数据
            break;
        }
        /*if (n % 10000 == 0)
        {
            ResultDeal();//处理数据
            PutOut();//输出数据
        }
        if (n >= 40000)
        {
            break;
        }*/
    }

    DeleteArray();//释放内存
    return 0;
}
// Implmentation of Liuping Wang MPC

#include "stdafx.h"
#include "iostream"
#include <windows.h>

#include "controller.h"

using namespace Eigen;
using namespace std;

int main()
{
	int N_sim = 100;
	MatrixXd xm(2, 1);
	xm << 0, 0;
	MatrixXd Xf = MatrixXd::Zero();

	system("pause");
    return 0;
}


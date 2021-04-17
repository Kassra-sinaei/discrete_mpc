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

	MPC my_mpc(10, 5, 0.1);
	// Plant State Space Model
	MatrixXd a(2, 2);
	MatrixXd b(2, 1);
	MatrixXd c(1, 2);
	a << 1,5,1,-5;
	b << 2, 1;
	c << 5, 7;
	my_mpc.setPlantModel(a, b, c);
	
	
	/*MatrixXd init_x(1, 1); // Plant initial state

	double ref = 10;

	for (int i = 0; i < N_sim; i++){

	}
	*/
	system("pause");
    return 0;
}


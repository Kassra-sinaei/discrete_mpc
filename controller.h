#pragma once

#include "stdafx.h"
#include "iostream"

#include"Eigen/Dense"
//#include "Eigen/eiquadprog.h"
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace std;

class MPC {
	friend int main();
	
public:
	MPC(int Np, int Nc);

private:
	MatrixXd model;
	int Np;
	int Nc;

	MatrixXd* c2d(MatrixXd A, MatrixXd B);
	MatrixXd* extrmodel(MatrixXd A, VectorXd B, VectorXd C);
	MatrixXd* mpcgains(MatrixXd Ap, MatrixXd Bp, MatrixXd Cp);
	double recede(MatrixXd Xf, MatrixXd F, MatrixXd phi, double rw, double rs);
};
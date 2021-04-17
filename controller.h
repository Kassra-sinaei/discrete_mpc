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
	MPC(int Np, int Nc, double timeStep);
	void setPlantModel(MatrixXd A, MatrixXd B, MatrixXd C);
private:
	double time_step;
	// MPC Controoler Gains for Prediction
	MatrixXd phi;
	MatrixXd F;

	// MPC Controoler Horizons
	int Np; // prediction Horizon
	int Nc; // control Horizon

	MatrixXd* c2d(MatrixXd A, MatrixXd B);
	MatrixXd* extrmodel(MatrixXd A, MatrixXd B, MatrixXd C);
	MatrixXd* mpcgains(MatrixXd Ap, MatrixXd Bp, MatrixXd Cp);
	double recede(MatrixXd Xf, MatrixXd F, MatrixXd phi, double rw, double rs);
};
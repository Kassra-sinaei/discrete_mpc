#include "stdafx.h"
#include "controller.h"

MPC::MPC(int Np, int Nc) {
	// Initialize controller parameters

	this->Nc = Nc;
	this->Np = Np;
}

MatrixXd* MPC::c2d(MatrixXd A, MatrixXd B) {
	MatrixXd ans[2];
	ans[0] = A.exp();
	const int dim = A.rows();
	ans[1] = (A.exp() - MatrixXd::Identity(dim, dim)) * B * A.transpose();
	return ans;
}

MatrixXd* MPC::extrmodel(MatrixXd A, VectorXd B, VectorXd C) {
	// Calculates augumented State Sapce matrices
	// Input Arguments are discretized state space matrices

	int m1 = C.rows();
	int n1 = C.cols();
	int n_in = B.cols();

	MatrixXd A_e = MatrixXd::Identity(n1 + m1, n1 + m1);
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n1; j++)
			A_e(i, j) = A(i, j);
	MatrixXd temp = C * A;
	for (int i = n1; i < n1 + m1; i++)
		for (int j = 0; j < n1; j++)
			A_e(i, j) = temp(i - n1, j);
	cout << "Ae: \n" << A_e << endl;

	MatrixXd B_e = MatrixXd::Zero(n1 + m1, n_in);
	temp = C * B;
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n_in; j++)
			B_e(i,j) = B(i,j);
	for (int i = n1; i < n1 + m1; i++)
		for (int j = 0; j < n_in; j++)
			B_e(i, j) = temp(i - n1, j);
	cout << "Be: \n" << B_e << endl;

	MatrixXd C_e = MatrixXd::Zero(m1, n1 + m1);
	temp = MatrixXd::Identity(m1, m1);
	for (int i = 0; i < m1; i++)
		for (int j = n1; j < n1 + m1; j++)
			C_e(i, j) = temp(i,j - n1);
	cout << "Ce: \n" << C_e << endl;

	MatrixXd ans[3];
	ans[0] = A_e;
	ans[1] = B_e;
	ans[2] = C_e;

	return ans;

}

MatrixXd* MPC::mpcgains(MatrixXd Ap, MatrixXd Bp, MatrixXd Cp) {
	// Calculates required components for prediction
	// Note that input Matrices should be augumented State Space matrices

	int m = Cp.cols();
	MatrixXd F(this->Np, m);

	VectorXd temp = Cp;
	for (int i = 0; i < this->Np; i++) {
		temp = temp * Ap;
		for (int j = 0; j < m; j++)
			F(i, j) = temp(0,j);
	}

	MatrixXd phi(this->Np, this->Nc);
	for (int i = 0; i < this->Np; i++) {
		for (int j = 0; j < this->Nc; j++) {
			if (i - j < 0)
				phi(i, j) = 0;
			if (i == j)
				phi(i, j) = (Cp * Bp)(0, 0);
			else
				phi(i, j) = (Cp * Ap.pow(i - j) * Bp)(0,0);
		}
	}

	MatrixXd ans[2];
	ans[0] = F;
	ans[1] = phi;
	return ans;
}

double MPC::recede(MatrixXd Xf, MatrixXd F, MatrixXd phi, double rw, double rs) {
	// Calculates Controller Output 
	// x(k), Y(k) --> delta_u
	// Xf: state feedback
	MatrixXd RS = rs * MatrixXd::Ones(this->Np,1);
	MatrixXd R_ = rw * MatrixXd::Identity(this->Nc, this->Nc);
	MatrixXd delta_u = (phi.transpose() * phi + R_).inverse() * phi.transpose() * (RS - F * Xf);
	return delta_u(0, 0);
}



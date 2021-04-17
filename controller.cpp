#include "stdafx.h"
#include "controller.h"

MPC::MPC(int Np, int Nc, double timeStep) {
	// Initialize controller parameters

	this->Nc = Nc;
	this->Np = Np;
	this->time_step = timeStep;
}

MatrixXd* MPC::c2d(MatrixXd A, MatrixXd B) {
	static MatrixXd ans[2];
	MatrixXd T = this->time_step * MatrixXd::Identity(A.rows(), A.cols());
	ans[0] = (A*T).exp();
	cout << "Discrete A:\n" << ans[0] << endl;
	const int dim = A.rows();
	ans[1] = A.inverse() * ((A*T).exp() - MatrixXd::Identity(dim, dim)) * B;
	cout << "Discrete B:\n" << ans[1] << endl;
	return ans;
}

void MPC::setPlantModel(MatrixXd A, MatrixXd B, MatrixXd C) {
	MatrixXd* discrete = (this->c2d(A, B));
	MatrixXd* augmented = this->extrmodel(discrete[0], discrete[1], C);
	MatrixXd* gains = this->mpcgains(augmented[0], augmented[1], augmented[2]);

	//this->phi = gains[1];
	//this->F = gains[0];
}

MatrixXd* MPC::extrmodel(MatrixXd A, MatrixXd B, MatrixXd C) {
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
	cout << "Augmented A: \n" << A_e << endl;

	MatrixXd B_e = MatrixXd::Zero(n1 + m1, n_in);
	temp = C * B;
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n_in; j++)
			B_e(i,j) = B(i,j);
	for (int i = n1; i < n1 + m1; i++)
		for (int j = 0; j < n_in; j++)
			B_e(i, j) = temp(i - n1, j);
	cout << "Augmented B: \n" << B_e << endl;

	MatrixXd C_e = MatrixXd::Zero(m1, n1 + m1);
	temp = MatrixXd::Identity(m1, m1);
	for (int i = 0; i < m1; i++)
		for (int j = n1; j < n1 + m1; j++)
			C_e(i, j) = temp(i,j - n1);
	cout << "Augmented C: \n" << C_e << endl;

	static MatrixXd ans[3];
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
	MatrixXd temp = Cp;
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
			else if (i == j)
				phi(i, j) = (Cp * Bp)(0, 0);
			else
				phi(i, j) = (Cp * Ap.pow(i - j) * Bp)(0,0);
		}
	}

	static MatrixXd ans[2];
	ans[0] = F;
	ans[1] = phi;
	cout << "Matrix F:\n" << ans[0] << endl;
	cout << "Matrix Phi:\n" << ans[1] << endl;
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

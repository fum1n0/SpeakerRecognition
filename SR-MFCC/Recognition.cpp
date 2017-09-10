#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include"Recognition.hpp"

Recognition::Recognition() { }

Recognition::Recognition(MFCC mfcc) {
	name = mfcc.name;
	ave = mfcc.ave;
	variance = mfcc.variance;
	determinant = mfcc.determinant;

	matrix = eigen_matrix(mfcc.variance); // C4996
	invertible = matrix.inverse();

	normarize = Recognition::calc_MLE(mfcc.ave);
	normarize = 1.0 / normarize;
	live = true;
}

bool Recognition::loadingModel() {

	if (!csv.isEmpty())csv.close();

	Optional<String> state;

	state = Dialog::GetOpenAll();
	if (!state.has_value()) {
		cout << "CSV Select Error" << endl << endl;
		return false;
	}
	csv.open(state.value());


	string check = csv.get<String>(0, 0).narrow();
	if (check != "8255850727") {
		cout << "Not Model CSV File" << endl << endl;
		return false;
	}

	name = csv.get<String>(1, 0).narrow();

	int size = csv.get<int32>(2, 0);

	ave = vector<double>(size, 0);
	variance = vector<vector<double>>(size, vector<double>(size, 0));

	for (int i = 0; i < size; i++) {
		ave[i] = csv.get<double>(3, i);
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			variance[i][j] = csv.get<double>(4 + i, j);
		}
	}

	determinant = csv.get<double>(14, 0);

	matrix = eigen_matrix(variance); // C4996
	invertible = matrix.inverse();

	normarize = Recognition::calc_MLE(ave);
	normarize = 1.0 / normarize;

	live = true;

	return true;
}


double Recognition::calc_MLE(vector<double>& X) {

	if (X.size() != ave.size())return -1;

	double theta = pow(2 * Pi, (double)ave.size() / 2.0);
	double r, tmp;

	vector<double> sub, tmp_vec;
	sub = X;
	tmp_vec = X;

	for (unsigned int i = 0; i < sub.size(); i++) {
		sub[i] -= ave[i];
	}

	for (unsigned int i = 0; i < ave.size(); i++) {
		tmp = 0;
		for (unsigned int j = 0; j < ave.size(); j++) {
			tmp += sub[i] * (double)invertible(i, j);
		}
		tmp_vec[i] = tmp;
	}

	r = 0;
	for (unsigned int i = 0; i < ave.size(); i++) {
		r += tmp_vec[i] * sub[i];
	}
	r *= -0.5;
	if (0 < r) {
		cout << "Mahalanobis' Distance Error" << endl;
		return -1;
	}

	r = exp(r) / (determinant*theta);

	if (isnan(r) || isinf(r)) {
		cout << "Calc Error" << endl;
		return -1;
	}

	if (&X == &ave) return r;
	else return normarize * r;
}

#pragma once


#include <Siv3D.hpp>
#include<vector>
#include<iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include<Eigen/LU>
#include "matrix_util.h"

#include"MFCC.hpp"

using namespace std;
using namespace Eigen;

class Recognition {

public:
	vector<double>ave;
	vector<vector<double>>variance;

	std::string name;
	Eigen::MatrixXf matrix;
	Eigen::MatrixXf invertible;
	double determinant;
	double normarize;
	CSVReader csv;

	bool live=false;

	Recognition();
	Recognition(MFCC);
	
	bool loadingModel();
	double calc_MLE(vector<double>&);

};
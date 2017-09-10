#pragma once

#include <Siv3D.hpp>
#include<vector>
#include<iostream>

using namespace std;

class Centroid {

public:
	
	vector<vector<double>>centroid;
	std::string name;
	double determinant;
	double normarize;
	CSVReader csv;
	bool live = false;


	bool loadingModel();
	double calc_MLE(vector<double>&);


};
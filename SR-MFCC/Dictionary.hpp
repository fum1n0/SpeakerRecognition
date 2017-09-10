#pragma once


#include <Siv3D.hpp>
#include<vector>
#include<iostream>
#include<fstream>
#include <functional>
#include <random>
#include"Centroid.hpp"

using namespace std;

class Dictionary {

public:
	vector<Centroid>dic;
	Centroid obs;
	vector<double>match;
	vector<double>euclidean;


	void load_dic(Centroid);
	void set_obs(Centroid);
	std::string calc_Recog();
	void dicReset();

};
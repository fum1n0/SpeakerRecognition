#pragma once

#include <Siv3D.hpp>
#include<vector>
#include<iostream>
#include<fstream>
#include <functional>
#include <random>

using namespace std;

class MFCC{

public:


	Sound sound;
	std::string name;

	Wave wav_base; // å≥îgå`
	Wave wav_pre; // çÇâπã≠í≤èàóùå„îgå`
	Wave wav_zero; // 0-pound

	unsigned long long int length;
	int Frame_L;
	int Frame_T;
	int leg;

	Wave wav_ana;

	int channels;

	std::vector<double>signal;

	std::vector<double>re;
	std::vector<double>im;
	std::vector<double>freq;

	vector<double>amp; // Amp
	vector<double>power; //power
	vector<double>dbv; // dBV
	
	vector<vector<double>>melFilterBank;
	vector<double>melFilterCenter;
	vector<double>melFilterStart;
	vector<double>melFilterEnd;
	vector<double>melFilterFreq;
	vector<double>ceps;

	vector<vector<double>>coefficients;

	int dimention = 12;
	vector<double>ave;
	vector<vector<double>>variance;
	double determinant;

	string fileMFCC;
	ofstream writing_MFCC;

	int num_codebook = 16;
	vector<int>cluster;
	vector<vector<double>>codebook;
	vector<vector<double>>delta;
	int epsilon = 2;
	vector<vector<double>>feature;

	/******************************************/

	MFCC(int, int);
	Wave erase_zeroAmp(Wave&);

	void init();

	void calc_MFCC(int);
	void hanning_execute(int);
	void fft_excute(vector<double>&, vector<double>&, vector<double>&, int);
	void calc_Normalization(vector<double>&);
	void calc_Power();
	void calc_Amp();

	double freqToMell(double);
	double melToFeeq(double);
	void calc_MelFilterBank();
	void MelFilter_execute();
	void DCT2_execute();

	
	void calc_VCM();
	void saveModelCSV();

	void kmeans(vector<vector<double>>);
	void calc_DeltaMFCC();
	void margeParameter();

};
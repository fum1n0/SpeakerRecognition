#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include"MFCC.hpp"


MFCC::MFCC(int frameL, int frameT) {

	while (sound.isEmpty() || sound.samplingRate() != 44100) {
		Optional<String> state = Dialog::GetOpenSound();
		if (!state.has_value()) {
			cout << "Sound Select Error" << endl << endl;
			return;
		}
		sound = Sound(state.value());
		Optional<AudioProperty> audio = Audio::GetProperty(state.value());
		name = audio->title.narrow();
		name = name.substr(0, name.size() - 4);
	}


	Frame_L = frameL;
	Frame_T = frameT;

	wav_zero = Wave(1); // 0-pound

	wav_base = sound.getWave();
	wav_base = erase_zeroAmp(wav_base);
	length = wav_base.lengthSample;
	leg = (int)(length / Frame_T);

	wav_zero[0] = Waving::DoubleToSample(0.0);

	for (int m = 0; m < Frame_L; m++)wav_base.push_back(wav_zero[0]);

	wav_pre = wav_base;
	wav_ana = Wave(Frame_L);

	// 高音強調
	for (int m = 1; m < length; m++)
		wav_pre[m] = Waving::DoubleToSample(((double)wav_base[m].left / 32768.0) - ((double)wav_base[m - 1].left / 32768.0) * 0.92);


	channels = 20;

	signal = vector<double>(Frame_L, 0);

	re = vector<double>(Frame_L, 0);;
	im = vector<double>(Frame_L, 0);
	freq = vector<double>(Frame_L, 0);;

	amp = vector<double>(Frame_L, 0); // Amp
	power = vector<double>(Frame_L, 0);//power
	dbv = vector<double>(Frame_L, 0); // dBV

	coefficients = vector<vector<double>>(dimention, vector<double>());
	delta = vector<vector<double>>(dimention, vector<double>());

	ave = vector<double>(dimention, 0);
	variance = vector<vector<double>>(dimention, vector<double>(dimention, 0.0));

	fileMFCC = name + "_MFCC.txt";;
	writing_MFCC.open(fileMFCC, ios::out);

	melFilterBank = vector<vector<double>>(channels, vector<double>(Frame_L / 2, 0));
	melFilterCenter = vector<double>(channels, 0);
	melFilterStart = vector<double>(channels, 0);
	melFilterEnd = vector<double>(channels, 0);
	melFilterFreq = vector<double>(channels, 0);


	MFCC::calc_MelFilterBank();


}


Wave MFCC::erase_zeroAmp(Wave& wav) {

	Wave erase;
	bool flag = true;
	for (size_t i = 1; i < wav.lengthSample - 1; i++) {
		flag = true;

		if (wav[i - 1].left == wav[i].left) {
			if (wav[i + 1].right == wav[i + 1].right) {
				flag = false;
			}
		}

		if (flag)erase.push_back(wav[i]);

	}

	erase.saveWAVE(L"erase.wav");

	return erase;


}




void MFCC::init() {
	signal = vector<double>(Frame_L, 0);
	re = vector<double>(Frame_L, 0);
	im = vector<double>(Frame_L, 0);
	freq = vector<double>(Frame_L, 0);
	amp = vector<double>(Frame_L, 0);
	power = vector<double>(Frame_L, 0);
	dbv = vector<double>(Frame_L, 0);

	melFilterFreq = vector<double>(channels, 0);
	ceps = vector<double>(channels, 0);
}




void MFCC::calc_MFCC(int indexFrame) {

	MFCC::init();
	MFCC::hanning_execute(indexFrame);
	MFCC::fft_excute(signal, re, im, 1);
	MFCC::calc_Amp();
	MFCC::calc_Power();
	MFCC::MelFilter_execute();
	MFCC::DCT2_execute();

}




void MFCC::hanning_execute(int indexFrame) {

	for (int k = 0; k < Frame_L; k++) {
		wav_ana[k] = Waving::DoubleToSample(((double)wav_pre[indexFrame*Frame_T + k].left / 32768.0) * (0.5 - 0.5*cos(2 * Pi*k / Frame_L)));
		signal[k] = wav_ana[k].left / 32768.0; // ハニング補正無し,正規化
	}
}


void MFCC::fft_excute(vector<double>& signal_ptr, vector<double>& re_ptr, vector<double>& im_ptr, int inv) {

	vector<double>sig = signal_ptr;

	double w;
	std::vector<double>wr;
	std::vector<double>wi;
	w = (double)inv * 2.0 * Pi / (double)sig.size();
	wr = std::vector<double>((int)sig.size(), 0);
	wi = std::vector<double>((int)sig.size(), 0);
	for (int i = 0; i < (int)sig.size(); i++) {
		wr[i] = cos(w * i);
		wi[i] = -sin(w*i);
	}


	long i = 0;

	if (inv > 0) {
		for (int j = 1; j < (int)sig.size() - 1; j++) {
			for (int k = (int)sig.size() >> 1; k > (i ^= k); k >>= 1); // k=2^{n-1}; k > i = i xor k; k/2
			if (j < i) {
				std::swap(sig[i], sig[j]);
			}
		}
		re_ptr = sig;

	}
	else {
		for (int j = 1; j < (int)sig.size() - 1; j++) {
			for (int k = (int)sig.size() >> 1; k > (i ^= k); k >>= 1); // k=2^{n-1}; k > i = i xor k; k/2
			if (j < i) {
				std::swap(re_ptr[i], re_ptr[j]);
				std::swap(im_ptr[i], im_ptr[j]);
			}
		}
	}


	double xr1, xr2, xi1, xi2;
	int s, t, turn1, turn2, turn3;

	for (int rep = 2; rep <= (int)sig.size(); rep *= 2) { // 一番外側のrep, 2DFT -> 4DFT -> ... -> NDFT
		s = (int)sig.size() / rep; // ブロック数
		t = rep / 2; // 1ブロックが保有する数の半分
		for (int u = 0; u < s; u++) { // NDFTの中のブロック数
			for (int v = 0; v < t; v++) { // ブロックの半分まで

				turn1 = v + rep*u;
				turn2 = s*v;
				turn3 = s*(v + t);

				xr1 = re_ptr[turn1];
				xr2 = re_ptr[turn1 + t];
				xi1 = im_ptr[turn1];
				xi2 = im_ptr[turn1 + t];


				re_ptr[turn1] = xr1 + xr2 * wr[turn2] - xi2 * wi[turn2];
				im_ptr[turn1] = xi1 + xi2 * wr[turn2] + xr2 * wi[turn2];

				re_ptr[turn1 + t] = xr1 + xr2 * wr[turn3] - xi2 * wi[turn3];
				im_ptr[turn1 + t] = xi1 + xi2 * wr[turn3] + xr2 * wi[turn3];


			}
		}
	}

	if (inv > 0) { // フーリエ変換
		for (int k = 0; k < (int)re_ptr.size(); k++) {
			//re_ptr[k] /= (double)re_ptr.size();
			//im_ptr[k] /= (double)im_ptr.size();
		}
	}
	else signal_ptr = re_ptr; // 逆フーリエ変換


}


void MFCC::calc_Normalization(vector<double>& num) {
	double max = *std::max_element(num.begin(), num.end());

	for (int i = 0; i < (int)num.size(); i++)num[i] /= max;
}


void MFCC::calc_Power() {

	for (int i = 0; i < amp.size(); i++) {
		power[i] = amp[i] * amp[i];
	}

}


void MFCC::calc_Amp() {

	for (int i = 0; i < amp.size(); i++) {
		amp[i] = sqrt(pow(re[i], 2.0) + pow(im[i], 2.0));
	}

}


double MFCC::freqToMell(double fr) {
	return 1127.01048 * log((fr / 700.0) + 1.0);
}


double MFCC::melToFeeq(double me) {
	return 700.0 * (exp(me / 1127.01048) - 1.0);
}


void MFCC::calc_MelFilterBank() {

	double melMax = MFCC::freqToMell((double)sound.samplingRate() / 2.0);
	double dmel = melMax / (double)(channels + 1);

	// フィルタバンク作成
	for (int i = 0; i < melFilterCenter.size(); i++)melFilterCenter[i] = (i + 1)*dmel;

	melFilterStart[0] = 0.0;
	for (int i = 1; i < melFilterStart.size(); i++)melFilterStart[i] = melFilterCenter[i - 1];

	melFilterEnd[melFilterEnd.size() - 1] = (melFilterEnd.size() + 1)*dmel;
	for (int i = 0; i < melFilterEnd.size() - 1; i++)melFilterEnd[i] = melFilterCenter[i + 1];

	for (int i = 0; i < melFilterCenter.size(); i++) { // メル→周波数→インデックス
		melFilterCenter[i] = (int)(Frame_L*melToFeeq(melFilterCenter[i]) / (double)sound.samplingRate());
		melFilterStart[i] = (int)(Frame_L*melToFeeq(melFilterStart[i]) / (double)sound.samplingRate());
		melFilterEnd[i] = (int)(Frame_L *melToFeeq(melFilterEnd[i]) / (double)sound.samplingRate());
	}

	for (int i = 0; i < melFilterBank.size(); i++) {
		double up = 1 / (double)(melFilterCenter[i] - melFilterStart[i]);
		double down = 1 / (melFilterEnd[i] - melFilterCenter[i]);

		for (int j = 0; j < melFilterBank[i].size(); j++) {
			if (melFilterStart[i] <= j && j <= melFilterCenter[i]) {
				melFilterBank[i][j] = (j - melFilterStart[i])*up;
			}
			else if (melFilterCenter[i] < j && j <= melFilterEnd[i]) {
				melFilterBank[i][j] = 1.0 - (j - melFilterCenter[i])*down;
			}
			else continue;
		}
	}


}


void MFCC::MelFilter_execute() {

	for (int i = 0; i < channels; i++) {

		for (int j = 0; j < melFilterBank[i].size(); j++) {
			melFilterFreq[i] += power[j] * melFilterBank[i][j];
		}
		if (melFilterFreq[i] != 0.0)melFilterFreq[i] = log10(melFilterFreq[i]);
		else melFilterFreq[i] = 0.0;
	}

}


void MFCC::DCT2_execute() {

	for (int j = 0; j < channels; j++)ceps[0] += melFilterFreq[j] / sqrt((double)channels);
	writing_MFCC << "0  " << ceps[0] << endl;

	//coefficients[0].push_back(ceps[0]);

	for (int i = 1; i < channels; i++) {
		for (int j = 0; j < channels; j++) {
			ceps[i] += melFilterFreq[j] * sqrt(2.0 / (double)channels) * cos(((2.0*j + 1.0)*i*Pi) / (2.0*(double)channels));
		}
		writing_MFCC << i << " " << ceps[i] << endl;
		if (4 <= i && i < 4 + dimention)coefficients[i - 4].push_back(ceps[i]);
		//if (i < dimention)coefficients[i].push_back(ceps[i]);
	}

}


void MFCC::calc_VCM() {
	int cnt;

	for (unsigned int s = 0; s < coefficients.size(); s++) {
		cnt = 0;

		for (int t = 0; t < (int)coefficients[s].size(); t++) {
			//if (coefficients[s][t] != 1) {
			cnt++;
			ave[s] += coefficients[s][t];
			//}
		}
		if (cnt != 0)ave[s] /= (double)cnt;
		else ave[s] = 1;

		for (int t = 0; t < (int)coefficients[s].size(); t++) {
			if (coefficients[s][t] != 1) {
				variance[s][s] += pow((ave[s] - coefficients[s][t]), 2);
			}
		}
		if (cnt != 0)variance[s][s] /= (double)cnt;
		else variance[s][s] = 0;

	}

	for (unsigned int i = 0; i < variance.size(); i++) {
		for (unsigned j = i; j < variance[i].size(); j++) {
			if (i == j)continue;
			cnt = 0;
			for (unsigned n = 0; n < coefficients[i].size(); n++) {
				if (coefficients[i][n] == 1 || coefficients[j][n] == 1)continue;
				variance[i][j] += ((ave[i] - coefficients[i][n])*(ave[j] - coefficients[j][n]));
				cnt++;
			}
			if (cnt != 0)variance[i][j] /= (double)cnt;
			else variance[i][j] = 0;
			variance[j][i] = variance[i][j];
		}
	}

	vector<vector<double>> copy;
	copy = variance;

	double det = 1.0, buf;

	//三角行列を作成
	for (unsigned int i = 0; i < copy.size(); i++) {
		for (unsigned int j = 0; j < copy[i].size(); j++) {
			if (i < j) {
				buf = copy[j][i] / copy[i][i];
				for (unsigned int k = 0; k < copy.size(); k++) {
					copy[j][k] -= copy[i][k] * buf;
				}
			}
		}
	}

	//対角部分の積
	for (unsigned int i = 0; i < copy.size(); i++) {
		det *= copy[i][i];
	}

	determinant = sqrt(fabs(det));

}


void MFCC::saveModelCSV() {

	writing_MFCC.close();

	string csv_name = "./csv/" + name + ".csv";

	CSVWriter writer(Widen(csv_name));
	writer.write(8255850727);
	writer.nextLine();
	writer.write(Widen(name));
	writer.nextLine();
	writer.write(ave.size());
	writer.nextLine();

	for (int i = 0; i < (int)ave.size(); i++)writer.write(ave[i]);
	writer.nextLine();

	for (int i = 0; i < (int)variance.size(); i++) {
		for (int j = 0; j < (int)variance[i].size(); j++) {
			writer.write(variance[i][j]);
		}
		writer.nextLine();
	}

	writer.write(determinant);

	cout << "Model Create Success: " << name << endl << endl;

}


void MFCC::kmeans(vector<vector<double>> sample) {

	if (sample.size() <= 0)return;

	cluster = vector<int>(sample[0].size());
	codebook = vector<vector<double>>(num_codebook, vector<double>(sample.size(), 0));
	vector<int>cnt(num_codebook, 0);

	vector<double>distance;
	vector<double>::iterator minPt;
	size_t minId;

	bool init = false;
	bool flag = false;
	random_device rnd;

	int num = 0;

	while (!flag) {
		flag = true;
		//num++;
		//cout << num << endl;
		// init
		if (init == false) {
			//cout << "init start" << endl;
			cnt = vector<int>(num_codebook, 0);
			codebook = vector<vector<double>>(num_codebook, vector<double>(sample.size(), 0));

			for (int i = 0; i < cluster.size(); i++) {
				cluster[i] = rnd() % num_codebook;
				for (int j = 0; j < sample.size(); j++)codebook[cluster[i]][j] += sample[j][i];
				cnt[cluster[i]]++;
			}

			for (int i = 0; i < codebook.size(); i++) {
				for (int j = 0; j < codebook[i].size(); j++)codebook[i][j] /= (double)cnt[i];
			}
			init = true;
		}


		// クラスタ更新
		for (int i = 0; i < sample[0].size(); i++) {
			distance = vector<double>(num_codebook, 0);
			for (int j = 0; j < codebook.size(); j++) {
				for (int k = 0; k < codebook[j].size(); k++) {
					distance[j] += pow(codebook[j][k] - sample[k][i], 2.0);
				}
				distance[j] = sqrt(distance[j]);
			}

			minPt = min_element(distance.begin(), distance.end());
			minId = std::distance(distance.begin(), minPt);

			if (cluster[i] != minId)flag = false;
			cluster[i] = minId;
		}

		// セントロイド更新
		cnt = vector<int>(num_codebook, 0);
		codebook = vector<vector<double>>(num_codebook, vector<double>(sample.size(), 0));

		for (int i = 0; i < cluster.size(); i++) {
			for (int j = 0; j < sample.size(); j++)codebook[cluster[i]][j] += sample[j][i];
			cnt[cluster[i]]++;
		}

		for (int i = 0; i < codebook.size(); i++) {
			if (cnt[i] == 0) {
				init = false;
				flag = false;
			}
			for (int j = 0; j < codebook[i].size(); j++)codebook[i][j] /= (double)cnt[i];
		}


	}


	// save
	string csv_name = "./csv/" + name + "_centroid.csv";
	CSVWriter writer(Widen(csv_name));
	writer.write(727825585);
	writer.nextLine();
	writer.write(Widen(name));
	writer.nextLine();
	writer.write(codebook[0].size());
	writer.nextLine();
	writer.write(num_codebook);
	writer.nextLine();

	for (int i = 0; i < (int)codebook.size(); i++) {
		for (int j = 0; j < (int)codebook[i].size(); j++) {
			writer.write(codebook[i][j]);
		}
		writer.nextLine();
	}

	cout << "Centroid Create Success: " << name << endl << endl;

}



void MFCC::calc_DeltaMFCC() {

	double denominator;
	double numerator = 0;

	for (int k = -1 * epsilon; k <= epsilon; k++)numerator += k*k;

	for (int i = 0; i < coefficients.size(); i++) {
		for (int j = 0; j < coefficients[i].size(); j++) {
			if (0 <= j - epsilon && j + epsilon < coefficients[i].size()) {
				denominator = 0.0;
				for (int k = -1 * epsilon; k <= epsilon; k++) {
					denominator += (double)k * coefficients[i][j + k];
				}

				delta[i].push_back(denominator / numerator);
			}
			else delta[i].push_back(0.0);
		}
	}


}

void MFCC::margeParameter() {

	string fileMarge;
	ofstream writing_Marge;
	fileMarge = name + "_marge.txt";;
	writing_Marge.open(fileMarge, ios::out);

	feature = vector<vector<double>>(coefficients.size() + delta.size(), vector<double>());

	for (int i = 0; i < coefficients[0].size(); i++) {

		if (0 <= i - epsilon && i + epsilon < coefficients[0].size()) {

			for (unsigned j = 0; j < coefficients.size(); j++) {
				feature[j].push_back(coefficients[j][i]);
				writing_Marge << j << " " << coefficients[j][i] << endl;
			}
			for (unsigned j = 0; j < delta.size(); j++) {
				feature[coefficients.size() + j].push_back(delta[j][i]);
				writing_Marge << coefficients.size() + j << " " << delta[j][i] << endl;
			}
		}


	}


	writing_Marge.close();

}
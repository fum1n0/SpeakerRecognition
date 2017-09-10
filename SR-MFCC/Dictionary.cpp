
#include"Dictionary.hpp"



void Dictionary::load_dic(Centroid cent) {
	dic.push_back(cent);
	cout << "Set Dictionary Model Name " << cent.name << endl;
}


void Dictionary::set_obs(Centroid cent) {
	obs = cent;
	cout << "Set Observer Model Name " << cent.name << endl;
}

std::string Dictionary::calc_Recog() {

	if (dic.size() <= 0 || obs.live == false)return "Not Set Centroid Model";

	
	match = vector<double>(dic.size(), 0);
	euclidean = vector<double>(obs.centroid.size(), 0);
	
	for (int i = 0; i < dic.size(); i++) {

		for (int j = 0; j < obs.centroid.size(); j++) {
			euclidean = vector<double>(obs.centroid.size(), 0);
			for (int k = 0; k < dic[i].centroid.size(); k++) {

				for (int l = 0; l < obs.centroid[j].size(); l++) {
					euclidean[k] += pow((obs.centroid[j][l] - dic[i].centroid[k][l]), 2.0);
				}

			}

			match[i] += *min_element(euclidean.begin(), euclidean.end());
		}

		cout << dic[i].name << " is match: " << match[i] << endl;
	}


	vector<double>::iterator minId = min_element(match.begin(), match.end());
	size_t minIndex = distance(match.begin(), minId);

	cout << "Result: "<< dic[minIndex].name << endl;
	


	return "Obs: " + obs.name + "\nResult: " + dic[minIndex].name;

}


void Dictionary::dicReset() {
	dic = vector<Centroid>();
	cout << "Reset Dictionary" << endl;
}

#include"Centroid.hpp"





bool Centroid::loadingModel() {

	if (!csv.isEmpty())csv.close();

	Optional<String> state;

	state = Dialog::GetOpenAll();
	if (!state.has_value()) {
		cout << "CSV Select Error" << endl << endl;
		return false;
	}
	csv.open(state.value());


	string check = csv.get<String>(0, 0).narrow();
	if (check != "727825585") {
		cout << "Not Centroid Model CSV File" << endl << endl;
		return false;
	}

	name = csv.get<String>(1, 0).narrow();

	int dimention = csv.get<int32>(2, 0);
	int num_code = csv.get<int32>(3, 0);

	centroid = vector<vector<double>>(num_code, vector<double>(dimention, 0));


	for (int i = 0; i < num_code; i++) {
		for (int j = 0; j < dimention; j++) {
			centroid[i][j] = csv.get<double>(4 + i, j);
		}
	}

	live = true;

	return true;
}
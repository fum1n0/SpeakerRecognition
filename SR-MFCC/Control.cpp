
#include "Control.hpp"


struct Matching : MyApp::Scene {

	double per = -2;
	String from_str;
	String to_str;
	String recog_name;

	void init() override {
		if (m_data->from.live && m_data->to.live)per = m_data->from.calc_MLE(m_data->to.ave);
		recog_name = Widen(m_data->dict.calc_Recog());
	}

	void update() override {

		if (!m_data->from.live) {
			from_str = L"From Model : Not Set";
			//FontAsset(L"Scene")(from_str).draw(10, 30);
		}
		else {
			from_str = L"From Model : ";
			//FontAsset(L"Scene")(from_str, Widen(m_data->from.name)).draw(10, 30);
		}


		if (!m_data->to.live) {
			to_str = L"To Model : Not Set";
			//FontAsset(L"Scene")(to_str).draw(10, 70);
		}
		else {
			to_str = L"To Model : ";
			//FontAsset(L"Scene")(to_str, Widen(m_data->to.name)).draw(10, 70);
		}

		if (per == -1) {
			cout << "Modelデータが正しくないため計算が行われませんでした" << endl << endl;
			per = -2;
		}
		else if (per != -2) {
			//FontAsset(L"Scene")(L"認識率 : ", DecimalPlace(2), per * 100, L"%").draw(10, 110);
		}

		int i=0;
		if (0 < m_data->dict.dic.size() && m_data->dict.obs.live) {
						
			int max_h = 0;
			for (int j = 0; j < m_data->dict.dic.size(); j++) {
				max_h = max(max_h ,FontAsset(L"Scene").region(Widen(m_data->dict.dic[j].name)).w);
			}


			for (i = 0; i < m_data->dict.dic.size(); i++) {
				FontAsset(L"Scene")(Widen(m_data->dict.dic[i].name)).draw(10, 30 * i);
				RectF(max_h + 20, 30 * i + 10, 400 - m_data->dict.match[i]*5, 20).draw(HSV(i));				
			}
		}
		FontAsset(L"Scene")(recog_name).draw(10, 30 * (i+1));


		if (Input::KeyC.clicked)changeScene(L"Clear", 0);
		if (Input::KeyF.clicked)changeScene(L"From Model Load", 0);
		if (Input::KeyT.clicked)changeScene(L"To Model Load", 0);
		if (Input::KeyM.clicked)changeScene(L"Model Create", 0);
		if (Input::KeyR.clicked)changeScene(L"Recording", 0);
		if (Input::KeyD.clicked)changeScene(L"DictionaryLoad", 0);
		if (Input::KeyO.clicked)changeScene(L"ObsSet", 0);

	}

};

struct ModelCreate : MyApp::Scene {

	MFCC* mfcc;

	int i = 0;
	void init() override {
		mfcc = new MFCC(m_data->Frma_L, m_data->Frame_T);

	}

	void update() override {

		if (Input::KeySpace.released) {
			cout << "モデル作成を中断しました" << endl << endl;
			changeScene(L"Matching", 0);
		}

		if (mfcc->sound.isEmpty()) { // 空の場合
			delete(mfcc);
			changeScene(L"Matching", 0);
			return;
		}

		mfcc->calc_MFCC(i);

		if (mfcc->leg <= i) {

			mfcc->calc_DeltaMFCC();
			mfcc->margeParameter();

			cout << "Clustering" << endl;
			mfcc->kmeans(mfcc->feature);

			mfcc->calc_VCM();
			mfcc->saveModelCSV();
			delete(mfcc);
			changeScene(L"Matching", 0);
		}

		//フレーム数表示
		FontAsset(L"Scene")(L"モデル作成: ", Widen(mfcc->name)).draw(10, 30);
		FontAsset(L"Scene")(L"進行フレーム数: ", i*m_data->Frame_T).draw(10, 70);
		int per = (int)((double)i*100.0 / (double)mfcc->leg);
		FontAsset(L"Scene")(L"解析進行度: ", per, L"%").draw(10, 110);

		i++;
	}

};

struct FromModelLoad : MyApp::Scene {
	void init() override {
		m_data->from.loadingModel();
	}

	void update() override { changeScene(L"Matching", 0); }

};

struct ToModelLoad : MyApp::Scene {
	void init() override {
		m_data->to.loadingModel();
	}

	void update() override { changeScene(L"Matching", 0); }

};

struct Recording : MyApp::Scene {

	int time = 10;
	string name, check;
	Recorder recorder;
	
	void init() override {

		cout << "10秒間録音をします" << endl;;
		
		cout << "保存ファイル名を入力してください: ";
		cin >> name;

		cout << "3秒後録音を開始します" << endl;
		System::Sleep(3500);

		if (!recorder.open(0, (int)(time * 44100))) {
			cout << "マイクデバイス エラー" << endl;
			return;
		}


	}

	void update() override {

		if (!recorder.start()) {
			cout << "録音開始 エラー" << endl;
			changeScene(L"Matching", 0);
			return;
		}
		cout << "録音開始" << endl;

		while (recorder.isEnded() != true) {}
		recorder.getWave().save(L"./record/" + Widen(name) + L".wav");
		cout << "録音終了" << endl << endl;
		changeScene(L"Matching", 0);

	}

};

struct DictionaryLoad : MyApp::Scene {
	void init() override {
		if (m_data->cent.loadingModel())m_data->dict.load_dic(m_data->cent);
	}
	void update() override { changeScene(L"Matching", 0); }
};


struct ObsSet : MyApp::Scene {
	void init() override {
		if (m_data->cent.loadingModel())m_data->dict.set_obs(m_data->cent);
	}
	void update() override { changeScene(L"Matching", 0); }
};

struct Clear : MyApp::Scene {
	void init() override {
		system("cls");
		m_data->dict.dicReset();

	}
	void update() override { changeScene(L"Matching", 0); }
};


Control::Control() {

	manager.setFadeColor(Palette::White);

	// シーンを設定
	manager.add<Matching>(L"Matching");
	manager.add<ModelCreate>(L"Model Create");
	manager.add<FromModelLoad>(L"From Model Load");
	manager.add<ToModelLoad>(L"To Model Load");
	manager.add<Recording>(L"Recording");
	manager.add<DictionaryLoad>(L"DictionaryLoad");
	manager.add<ObsSet>(L"ObsSet");
	manager.add<Clear>(L"Clear");
}
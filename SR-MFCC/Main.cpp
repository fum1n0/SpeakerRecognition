#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include<Siv3D.hpp>
#include <HamFramework.hpp>

#include "Control.hpp"

using namespace std;
using namespace Eigen;


void Main() {
	ScalableWindow::Setup(600, 480);
	Window::SetTitle(L"Speaker Recognition");
	FontAsset::Register(L"Scene", 20);
	Control cont;

	Console::Open(); // Console Window
		
	while (System::Update()) {
		if (!cont.manager.updateAndDraw())break;		
	}

	Console::Close();
}
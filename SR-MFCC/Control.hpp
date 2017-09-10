#pragma once
#include<Siv3D.hpp>
# include <HamFramework.hpp>

#include<iostream>
#include <limits>

#include "MFCC.hpp"
#include "Recognition.hpp"
#include "Centroid.hpp"
#include "Dictionary.hpp"

struct ModelData{
	Recognition from;
	Recognition to;
	Centroid cent;
	Dictionary dict;

	int Frma_L = 4096;
	int Frame_T = 256;

};

using MyApp = SceneManager<String, ModelData>;


class Control {

public:
	
	MyApp manager;

	Control();
	
};
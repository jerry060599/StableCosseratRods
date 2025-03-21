#pragma once
#include "Cosserat.h"

namespace Cosserat {
	Sim* generateCantileverExample();
	Sim* generateVBDFailureExample();
	Sim* generateTreeExample(int numTrees = 32);
	Sim* generateLargeBridgeExample();
	Sim* generateSmallBridgeExample();
	Sim* generateSlingshotExample();
	Sim* generateSlinkyExample();
}
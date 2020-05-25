#pragma once
#include "BasisSet.h"

class HFCalculator
{
public:

	void init();

	BasisSet& getSet();

private:

	BasisSet _bset;


};
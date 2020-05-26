#pragma once
#include "Mol.h"
#include <string>

class HFCalculator
{
public:

	static void run();
	static std::string nucleonsNames[10];
	static double bohrAngstrom;

private:

	static int atomFromConsole();


};
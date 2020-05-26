#include "HFCalculator.h"

std::string HFCalculator::nucleonsNames[10] = {
	"H", "He", "Li", "Be", "B", "C", "N","O","F","Ne"
};

double HFCalculator::bohrAngstrom = 0.52917721;

void HFCalculator::run()
{
	int x = 1;
	while (x)
	{
		std::cout << "Let's start some quantum chemistry computations..." << std::endl;
		std::cout << "Give number of atoms in molecule: ";
		int n;
		std::cin >> n;
		std::cout << std::endl;
		std::cout << "Atoms to choose from: [ ";
		for (int i = 0; i < 9; i++)
		{
			std::cout << nucleonsNames[i] << ", ";
		}
		std::cout << nucleonsNames[9] << " ]"<<std::endl;
		Mol molecule;
		for (int i = 0; i < n; i++)
		{
			std::cout << "Choose " << i + 1 << " atom: ";
			int e;
			Position p;
			e = atomFromConsole();
			while (e == -1)
			{
				std::cout << "Choose one atom from list" << std::endl;
				e = atomFromConsole();
			}
			std::cout << "Give atom " << i + 1 << " position in angstroms: ";
			std::cin >> p.x;
			std::cin >> p.y;
			std::cin >> p.z;
			p.x /= bohrAngstrom;
			p.y /= bohrAngstrom;
			p.z /= bohrAngstrom;
			p.reloadSpherical();
			molecule.addNucleon(Nucleon{ p, e });

		}
		std::cout << "Molecule created" << std::endl;
		std::cout << "Loading basis set..." << std::endl;
		molecule.initBasisSet();
		std::cout << "Calculating integrals..." << std::endl;
		molecule.calculateIntegrals();
		std::cout << "HF Self Consistent Procedure..." << std::endl;
		molecule.HF_TO_Divergance(0.001);
		std::cout << "And here we have the result" << std::endl;
		std::cout << "Molecule energy is equal: " << molecule.getMoleculeEnergy() << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "Once again? 1-yes/0-no" << std::endl;
		std::cin >> x;
	}
}

int HFCalculator::atomFromConsole()
{
	std::string name;
	std::cin >> name;
	for (size_t i = 0; i < 10; i++)
	{
		if (name == nucleonsNames[i])
			return i + 1;
	}
	return -1;
}

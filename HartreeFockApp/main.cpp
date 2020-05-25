#include <iostream>

#include "Mol.h"

int main()
{

	Mol O{ 8 };
	Mol C{ 6 };
	O.initBasisSet();
	O.calculateIntegrals();
	O.HF_TO_Divergance(0.001);
	C.initBasisSet();
	C.calculateIntegrals();
	C.HF_TO_Divergance(0.001);

	Mol H20{ { {{0,0,0}, 8}, {{-1.8897,0,0 }, 1}, {{1.8897 * 0.2588190451,1.8897 * 0.96592582628906 ,0  },1} } };
	H20.initBasisSet();
	H20.calculateIntegrals();
	H20.HF_TO_Divergance(0.001);
	std::cout << H20.getMoleculeEnergy() << std::endl;
	std::cout << "Jonization energy:  " << O.getMoleculeEnergy() - 1.0 - H20.getMoleculeEnergy() << std::endl;

	Mol CO2{{ {{-2.19775,0,0 }, 8} , {{0,0,0}, 6}, {{2.19775,0,0 }, 8} } };

	CO2.initBasisSet();
	CO2.calculateIntegrals();
	CO2.HF_TO_Divergance(0.001);
	std::cout << CO2.getMoleculeEnergy() << std::endl;
	std::cout << "Jonization energy:  " << 2*O.getMoleculeEnergy() +C.getMoleculeEnergy() - CO2.getMoleculeEnergy() << std::endl;

	return 0;
}
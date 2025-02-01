#include <iostream>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include "madsen_root_finder.hpp"

int main()
{
	std::vector<double> polynomial{ 1, 2, 3, 4, 5, 6, 7, 8 };

	bitlush::madsen_root_finder<double> madsen;

	madsen.find_roots(polynomial.data(), polynomial.size() - 1);

	for (std::size_t i = 0; i < polynomial.size() - 1; i++)
	{
		std::cout << std::setprecision(15) << "root: " << madsen[i] << std::endl;
	}

	return EXIT_SUCCESS;
}
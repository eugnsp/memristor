#include "simulator.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>

int main(int argc, const char** argv)
{

	es_fe::Point pt1{4, 5};

	std::cout << pt1;

	try
	{
		std::cout << "Memristor simulation code ver. 1.0\n"
				  << "==================================\n"
				  << std::endl;

		Simulator sim;
		sim.run(argc, argv);
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception!\n" << e.what() << '\n';
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cout << "Exception!\n";
		return EXIT_FAILURE;
	}

	std::cout << "Simulation completed.\n";
	return EXIT_SUCCESS;
}

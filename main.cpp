#include <iostream>
#include <fstream>
#include <string>
#include "psstFiniteVolumeGrid.h"
#include "psstGridReader.h"
#include "psstAnsysGridReader.h"
#include "psstOutputManager.h"

int main(int argc, char*argv[])
{
	//problem_manager reading;
	const std::string problem_manager("problem_manager");
	std::ifstream in_file(problem_manager);

	std::cout << "Opening problem_manager file..." << std::endl;
	if (!(in_file.is_open()))
	{
		std::cerr << "Error: Could not open problem_manager file!";
		return 1;
	}
	std::string grid_name;
	std::string output_name;
	std::cout << "Reading grid name..." << std::endl;
	in_file >> grid_name;
	std::cout << "Grid name: " << grid_name << std::endl;

	std::cout << "Reading output files name..." << std::endl;
	in_file >> output_name;
	std::cout << "Name for output files: " << output_name;
	in_file.close();

	std::cout << "Problem_manager has been successfully read!" << std::endl;

	//create input grid;
	in_file.open(grid_name);
	if (!(in_file.is_open()))
	{
		std::cerr << "Error: Could not open \"" << grid_name << "\" file!";
		return 1;
	}

	std::cout << "\nReading grid..." << std::endl;
	psstInterfaceGridBuilder * builder = new psstAnsysBuilder;
	psstInterfaceGridDirector * director = new psstAnsysDirector;
	director->set_builder(builder);
	psstInterfaceInputGrid * grid = director->construct_grid(in_file);
	in_file.close();

	std::shared_ptr<psstInterfaceOutputManager> output_manager = std::make_shared<psstAnsysOutputManager>(output_name);
	psstGridFV * finite_volume_grid = grid->create_finite_volume_grid();
	finite_volume_grid->attach_output_manager(output_manager);
	finite_volume_grid->settings_file_name(output_name);
	finite_volume_grid->solver_parameters();
	finite_volume_grid->explicit_solve();

	delete finite_volume_grid;
	delete grid;
	delete director;
	delete builder;
	
	std::cout << "Done!" << std::endl;
	return 0;
}
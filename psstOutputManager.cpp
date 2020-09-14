#include <fstream>
#include <iomanip>
#include "psstOutputManager.h"
#include "psstFiniteVolumeGrid.h"
#include "psstFiniteVolumeGridTools.h"
#include "psstVariable.h"

void psstAnsysOutputManager::initialize()
{
	bc_types_names_[1] = "fluid";
	bc_types_names_[2] = "interior";
	bc_types_names_[3] = "wall";
	bc_types_names_[4] = "pressure-inlet";
	bc_types_names_[5] = "pressure-outlet";
	bc_types_names_[7] = "symmetry";
	bc_types_names_[8] = "periodic-shadow";
	bc_types_names_[9] = "pressure-far-field";
	bc_types_names_[10] = "velocity-inlet";
	bc_types_names_[12] = "periodic";
	bc_types_names_[14] = "radiator";
	bc_types_names_[20] = "mass-flow-inlet";
	bc_types_names_[24] = "interface";
	bc_types_names_[31] = "parent";//face
	bc_types_names_[32] = "parent";//cell
	bc_types_names_[36] = "outflow";
	bc_types_names_[37] = "axis";
}

psstInterfaceOutputManager * psstAnsysOutputManager::clone() const
{
	psstInterfaceOutputManager * new_manager = new psstAnsysOutputManager(*this);
	return new_manager;
}

void psstAnsysOutputManager::output_data(const psstGridFV & grid) const
{
	private_write_field(grid);
	private_write_cas_file(grid);
}

double psstAnsysOutputManager::private_get_variable_by_id(std::shared_ptr<const psstInterfaceVolume> vol, size_t variable_id) const
{
	double variable = 0.0;
	switch (variable_id)
	{
	case pressure_id:
		variable = vol->variable()->pressure();
		break;
	case density_id:
		variable = vol->variable()->density();
		break;
	case temperature_id:
		variable = vol->variable()->temperature();
		break;
	case vel_x_id:
		variable = vol->variable()->velocity(0);
		break;
	case vel_y_id:
		variable = vol->variable()->velocity(1);
		break;
	case vel_z_id:
		variable = vol->variable()->velocity(2);
		break;
	default:
		throw std::invalid_argument("Error: Unknown variable id!");
		break;
	}
	return variable;
}

void psstAnsysOutputManager::private_write_variable_to_tecplot(const std::map<size_t, psstVolumeStorage> * storage, size_t variable_id) const
{
	std::string file_name = file_name_ + ".dat";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);
	out_file.precision(11);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	size_t first_index = 0;
	size_t last_index = 0;
	double variable = 0.0;
	std::map<size_t, psstVolumeStorage>::const_iterator iter = storage->cbegin();
	for (iter; iter != storage->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ != 32)
		{
			first_index = last_index;
			last_index += iter->second.volumes_.size();
			out_file << "(\t " << 300 << " \t(\t " << variable_id << " \t " << iter->first
				<< " \t " << 1 << " \t " << 0 << " \t " << 0 << " \t " << first_index + 1 << " \t " << last_index << " \t) "
				<< std::endl << "(\n";

			for (size_t j = 0; j < iter->second.volumes_.size(); ++j)
			{
				variable = private_get_variable_by_id(iter->second.volumes_[j], variable_id);
				out_file << std::fixed << variable << "\n";
			}

			out_file << ")\n" << std::endl;
		}
	}

	out_file.close();
}

void psstAnsysOutputManager::private_write_field(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".dat";
	std::ofstream out_file;
	out_file.open(file_name);
	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");
	out_file.close();

	const std::map<size_t, psstVolumeStorage> * storage = grid.volume_storage();
	
	private_write_variable_to_tecplot(storage, pressure_id);
	private_write_variable_to_tecplot(storage, temperature_id);
	private_write_variable_to_tecplot(storage, density_id);
	private_write_variable_to_tecplot(storage, vel_x_id);
	private_write_variable_to_tecplot(storage, vel_y_id);
	if (grid.dimension() == 3)
	{
		private_write_variable_to_tecplot(storage, vel_z_id);
	}
}

void psstAnsysOutputManager::private_write_dimension(size_t dim) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	out_file << "(0 \"Dimension:\")\n";
	out_file << "(" << write_dimension << " " << dim << ")\n" << std::endl;

	out_file.close();
}

void psstAnsysOutputManager::private_write_vertices(const std::vector<std::shared_ptr<psstInterfaceVertex>> * vertices) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");


	size_t zone_id = 0;
	size_t first_index = 1;
	size_t last_index = vertices->size();
	size_t type = 1;
	size_t dimension = vertices->at(0)->dimension();

	out_file << "(0 \"Vertices:\")\n";
	out_file << "(" << std::dec << write_vertices << " (" << zone_id << " " << std::hex << first_index << " " << last_index << " " << type << " " << dimension << "))\n";
	zone_id = 1;
	out_file << "(" << std::dec << write_vertices << " (" << zone_id << " " << std::hex << first_index << " " << last_index << " " << type << " " << dimension << ")(\n";

	for (size_t i = 0; i < vertices->size(); ++i)
	{
		for (size_t j = 0; j < dimension; ++j)
			out_file << std::dec << std::setprecision(20) << std::fixed << std::left << std::scientific << " " << vertices->at(i)->component(j) << " ";
		out_file << "\n";
	}

	out_file << "))\n" << std::endl;
	out_file.close();
}

void psstAnsysOutputManager::private_write_faces(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	out_file << "(0 \"Faces:\")\n";

	size_t zone_id = 0; 
	size_t temp = 1;
	out_file << "(" << std::dec << write_faces << " (" << std::hex << zone_id << " " << temp << " " << grid.number_of_faces() << " " << 0 << "))\n";

	const std::map<size_t, psstFaceStorage> * face_storage = grid.face_storage();
	std::map<size_t, psstFaceStorage>::const_iterator iter = face_storage->cbegin();
	for (iter; iter != face_storage->cend(); ++iter)
	{
		
		zone_id = iter->first;
		out_file << "(" << std::dec << write_faces << " (" << std::hex << iter->first << " " << iter->second.first_index_ << " "
			<< iter->second.last_index_ << " " << iter->second.bc_type_id_ << " " << 0 << ")(\n";
		
		for (size_t i = 0; i < iter->second.faces_.size(); ++i)
		{
			out_file << std::hex << std::left << " " << iter->second.faces_.at(i)->type() << " ";
			for (size_t j = 0; j < iter->second.faces_.at(i)->number_of_vertices(); ++j)
				out_file << std::hex << std::left << " " << iter->second.faces_.at(i)->vertex(j)->id() << " ";
			for (size_t j = 0; j < 2; ++j)
			{
				if (iter->second.faces_.at(i)->neighbour(j).get() == 0)
				{
					out_file << std::hex << std::left << " " << 0 << " ";
				}
				else
				{
					out_file << std::hex << std::left << " " << iter->second.faces_.at(i)->neighbour(j)->id() << " ";
				}
			}
			out_file << std::endl;
		}
	
		out_file << "))\n" << std::endl;
	}

	out_file.close();
}

void psstAnsysOutputManager::private_write_cells(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	out_file << "(0 \"Cells:\")\n";

	size_t zone_id = 0;
	out_file << "(" << std::dec << write_cells << " (" << std::hex << zone_id << " " << 1 << " " << grid.number_of_inner_volumes() << " " << 0 << "))\n";

	const std::map<size_t, psstVolumeStorage> * volume_storage = grid.volume_storage();
	std::map<size_t, psstVolumeStorage>::const_iterator iter = volume_storage->cbegin();
	for (iter; iter != volume_storage->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 32 || iter->second.bc_type_id_ == 1)
		{

			zone_id = iter->first;
			out_file << "(" << std::dec << write_cells << " (" << std::hex << zone_id << " " << iter->second.first_index_ << " " << iter->second.last_index_ << " " << iter->second.bc_type_id_ << " " << 0 << ")(\n";
			for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
			{
				if (i % 9 == 0 && i != 0)
					out_file << std::endl;
				out_file << " " << iter->second.volumes_.at(i)->type();
			}
			out_file << std::endl << "))\n" << std::endl;
		}
	}

	out_file.close();
}

void psstAnsysOutputManager::private_write_cell_tree(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	const std::map<size_t, psstVolumeStorage> * volume_storage = grid.volume_storage();
	std::map<size_t, psstVolumeStorage>::const_iterator iter = volume_storage->cbegin();

	for (iter; iter != volume_storage->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ != 32)
			continue;
		out_file << "(0 \"Cells Trees:\")\n";
		break;
	}

	for (iter; iter != volume_storage->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ != 32)
			continue;
		out_file << "(" << std::dec << write_cell_tree << " (" << std::hex << iter->second.first_index_ << " " << iter->second.last_index_ << " "
			<< iter->first << " " << iter->second.children_zone_id_ << ")(\n";
		
		for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
		{
			out_file << " " << iter->second.volumes_.at(i)->number_of_children();
			for (size_t j = 0; j < iter->second.volumes_.at(i)->number_of_children(); ++j)
			{
				out_file << " " << iter->second.volumes_.at(i)->child(j)->id();
			}
			out_file << "\n";
		}

		out_file << "))\n" << std::endl;
	}

	out_file.close();
}

void psstAnsysOutputManager::private_write_face_tree(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	const std::map<size_t, psstFaceStorage> * face_storage = grid.face_storage();
	std::map<size_t, psstFaceStorage>::const_iterator iter = face_storage->cbegin();

	for (iter; iter != face_storage->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ != 31)
			continue;
		out_file << "(0 \"Faces Trees:\")\n";
		break;
	}

	for (iter; iter != face_storage->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ != 31)
			continue;
		out_file << "(" << std::dec << write_face_tree << " (" << std::hex << iter->second.first_index_ << " " << iter->second.last_index_ << " "
			<< iter->first << " " << iter->second.children_zone_id_ << ")(\n";

		for (size_t i = 0; i < iter->second.faces_.size(); ++i)
		{
			out_file << " " << iter->second.faces_.at(i)->number_of_children();
			for (size_t j = 0; j < iter->second.faces_.at(i)->number_of_children(); ++j)
			{
				out_file << " " << iter->second.faces_.at(i)->child(j)->id();
			}
			out_file << "\n";
		}

		out_file << "))\n" << std::endl;
	}

	out_file.close();
}

void psstAnsysOutputManager::private_write_zones(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name, std::ios::app);

	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.dat file!");

	out_file << "(0 \"Zones:\")\n";
	const std::map<size_t, psstVolumeStorage> * volume_storage = grid.volume_storage();
	std::map<size_t, psstVolumeStorage>::const_iterator vol_iter = volume_storage->cbegin();
	for (vol_iter; vol_iter != volume_storage->cend(); ++vol_iter)
	{
		if (vol_iter->second.bc_type_id_ != 1)
			continue;
		out_file << "(" << std::dec << write_zone_39 << " (" << vol_iter->first << " " << bc_types_names_.at(vol_iter->second.bc_type_id_) << " " << vol_iter->second.name_ << ")())\n" << std::endl;
	}

	const std::map<size_t, psstFaceStorage> * face_storage = grid.face_storage();
	std::map<size_t, psstFaceStorage>::const_iterator face_iter = face_storage->cbegin();
	for (face_iter; face_iter != face_storage->cend(); ++face_iter)
	{
		if (face_iter->second.bc_type_id_ == 31)
			continue;
		out_file << "(" << std::dec << write_zone_39 << " (" << face_iter->first << " " << bc_types_names_.at(face_iter->second.bc_type_id_) << " " << face_iter->second.name_ << ")())\n" << std::endl;
	}

	out_file.close();
}

void psstAnsysOutputManager::private_write_cas_file(const psstGridFV & grid) const
{
	std::string file_name = file_name_ + ".cas";
	std::ofstream out_file;
	out_file.open(file_name);
	if (!(out_file.is_open()))
		throw std::exception("Error: Could not open output *.cas file!");
	out_file.close();

	size_t dimension = grid.dimension();
	const std::vector<std::shared_ptr<psstInterfaceVertex>> * vertices = grid.vertex_storage();

	private_write_dimension(dimension);
	private_write_vertices(vertices);
	private_write_faces(grid);
	private_write_cells(grid);
	private_write_cell_tree(grid);
	private_write_face_tree(grid);
	private_write_zones(grid);
}
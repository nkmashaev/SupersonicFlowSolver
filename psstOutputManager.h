#include <string>
#include <map>
#include <memory>
#include <vector>
#ifndef _PSST_OUTPUT_MANAGER_H_
#define _PSST_OUTPUT_MANAGER_H_

class psstGridFV;
class psstInterfaceVolume;
struct psstVolumeStorage;
class psstInterfaceVertex;
class psstInterfaceOutputManager
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceOutputManager * clone() const = 0;
	virtual void output_data(const psstGridFV &) const = 0;
	virtual ~psstInterfaceOutputManager() {};
};

class psstInterfaceVolume;
class psstAnsysOutputManager : public psstInterfaceOutputManager
{
private:
	std::string file_name_;
	std::map<size_t, std::string> bc_types_names_;

	enum ansys_variable_data_id
	{
		pressure_id = 1,
		temperature_id = 3,
		density_id = 101,
		vel_x_id = 111,
		vel_y_id = 112,
		vel_z_id = 113
	};

	enum ansys_operation_id
	{
		write_dimension = 2,
		write_vertices = 10,
		write_cells = 12,
		write_faces = 13,
		write_zone_39 = 39,
		write_zone_45 = 45,
		write_cell_tree = 58,
		write_face_tree = 59
	};

	double private_get_variable_by_id(std::shared_ptr<const psstInterfaceVolume> vol, size_t variable_id) const;
	void private_write_variable_to_tecplot(const std::map<size_t, psstVolumeStorage> * storage, size_t variable_id) const;
	void private_write_field(const psstGridFV & grid) const;
	
	void private_write_dimension(size_t dimension) const;
	void private_write_vertices(const std::vector<std::shared_ptr<psstInterfaceVertex>> * vertices) const;
	void private_write_faces(const psstGridFV & grid) const;
	void private_write_cells(const psstGridFV & grid) const;
	void private_write_cell_tree(const psstGridFV & grid) const;
	void private_write_face_tree(const psstGridFV & grid) const;
	void private_write_zones(const psstGridFV & grid) const;
	void private_write_cas_file(const psstGridFV & grid) const;
public:
	psstAnsysOutputManager() :
		file_name_()
	{
		initialize();
	}

	psstAnsysOutputManager(const std::string file_name) :
		file_name_(file_name)
	{
		initialize();
	};

	virtual void initialize();
 	virtual psstInterfaceOutputManager * clone() const;

	virtual void output_data(const psstGridFV &) const;
	virtual ~psstAnsysOutputManager() {};
};

#endif
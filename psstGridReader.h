#ifndef _PSST_GRID_READER_H_
#define _PSST_GRID_READER_H_

#include<vector>
#include "psstFiniteVolumeGrid.h"

class psstInterfaceSimpleVertex
{
public:
	virtual psstInterfaceSimpleVertex * clone() const = 0;
	virtual ~psstInterfaceSimpleVertex() {};
};

class psstInterfaceSimpleFace
{
public:
	virtual psstInterfaceSimpleFace * clone() const = 0;
	virtual ~psstInterfaceSimpleFace() {};
};

class psstInterfaceSimpleCell
{
public:
	virtual psstInterfaceSimpleCell * clone() const = 0;
	virtual ~psstInterfaceSimpleCell() {};
};

class psstInterfaceSimpleStorage
{
public:
	virtual psstInterfaceSimpleStorage * clone() const = 0;
	virtual ~psstInterfaceSimpleStorage() {};
};

class psstInterfaceSimpleConnection
{
public:
	virtual psstInterfaceSimpleConnection * clone() const = 0;
	virtual ~psstInterfaceSimpleConnection() {};
};

class psstInterfaceInputGrid
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceInputGrid * clone() const = 0;
	virtual size_t dimension() const = 0;
	virtual void set_dimension(size_t i) = 0;

	virtual void reserve_vertices(size_t i) = 0;
	virtual size_t vertices_number() const = 0;
	virtual void add_vertex(const double *, size_t) = 0;
	virtual void add_vertex(const std::vector<double> &) = 0;
	virtual void add_vertex(const psstInterfaceSimpleVertex & vertex) = 0;

	virtual void reserve_faces(size_t face_numb) = 0;
	virtual size_t faces_number() const = 0;
	virtual void add_face(const size_t * vertex_index, size_t vol_left, size_t vol_right, size_t type, size_t zone_id) = 0;
	virtual void add_face(const std::vector<size_t> & vertex_index, size_t vol_left, size_t vol_right, size_t type, size_t zone_id) = 0;
	virtual void add_face(const psstInterfaceSimpleFace & face) = 0;
	virtual void add_face_connection(size_t * children, size_t number_of_children, size_t parent_index) = 0;

	virtual void reserve_cells(size_t face_numb) = 0;
	virtual size_t cells_number() const = 0;
	virtual void add_cell(size_t type, size_t zone_id) = 0;
	virtual void add_cell(const psstInterfaceSimpleCell & cell) = 0;
	virtual void add_cell_connection(size_t * children, size_t number_of_children, size_t parent_index) = 0;

	virtual size_t number_of_zones() const = 0;
	virtual void add_zones_storage(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type) = 0;
	virtual void add_zones_storage(size_t zone_id, const std::string & name) = 0;
	virtual void add_zones_storage(const psstInterfaceSimpleStorage & zone_storage) = 0;
	virtual void add_zones_storage_connection(size_t first_index, size_t last_index, size_t parent_zone_id, size_t children_zone_id) = 0;

	virtual psstGridFV * create_finite_volume_grid() const = 0;

	virtual ~psstInterfaceInputGrid() {};
};

class psstInterfaceGridBuilder
{
public:

	virtual void reset() = 0;
	virtual psstInterfaceInputGrid * get_grid() const = 0;
	virtual void read_dimension(size_t dim) const = 0;
	virtual void read_vertices(size_t zone_id, size_t first_index, size_t last_index, size_t dimension, std::istream * input_stream) const = 0;
	virtual void read_faces(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type, size_t face_type, std::istream * input_stream) const = 0;
	virtual void read_cells(size_t zone_id, size_t fisrt_index, size_t last_index, size_t bc_type, size_t cell_type, std::istream * input_stream) const = 0;
	virtual void read_zones(size_t zone_id, const std::string zone_name, std::istream * input_stream) const = 0;
	virtual void read_face_tree(size_t first_index, size_t last_index, size_t parent_id, size_t children_id, std::istream * input_stream) const = 0;
	virtual void read_cell_tree(size_t first_index, size_t last_index, size_t parent_id, size_t children_id, std::istream * input_stream) const = 0;
	virtual ~psstInterfaceGridBuilder() {};
};

class psstInterfaceGridDirector
{
public:
	virtual void set_builder(psstInterfaceGridBuilder * builder) = 0;
	virtual size_t operation_id() const = 0;
	virtual void set_operation_id(size_t operation_id) = 0;
	virtual void build_step() = 0;
	virtual psstInterfaceInputGrid * construct_grid(std::ifstream &) = 0;
	virtual ~psstInterfaceGridDirector() {};
};
#endif

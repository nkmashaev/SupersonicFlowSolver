#ifndef _PSST_ANSYS_GRID_READER_H_
#define _PSST_ANSYS_GRID_READER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "psstGridReader.h"
#include "psstFiniteVolumeGridTools.h"

class psstSimpleVertex : public psstInterfaceSimpleVertex
{
private:
	double * x_;
	size_t dimension_;

	psstSimpleVertex() {};
	void private_set_dimension(size_t dim);
public:
	psstSimpleVertex(size_t dim);
	psstSimpleVertex(const psstSimpleVertex &);
	psstSimpleVertex(const double *, size_t);
	psstSimpleVertex(const std::vector<double> &);

	virtual void initialize();
	virtual psstInterfaceSimpleVertex * clone() const;
	virtual void copy(const psstSimpleVertex &);

	virtual double * components() { return x_; };
	virtual const double * components() const { return x_; };
	virtual size_t dimension() const { return dimension_; };

	virtual double & component(size_t i) { return const_cast<double &>(static_cast<const psstSimpleVertex &>(*this).component(i)); };
	virtual const double & component(size_t i) const;

	psstSimpleVertex & operator=(const psstSimpleVertex &);
	virtual ~psstSimpleVertex() { delete[]x_; };
};

class psstSimpleFace : public psstInterfaceSimpleFace
{
private:
	size_t zone_id_;
	size_t type_;
	size_t vertex_number_;
	size_t * vertex_index_;
	size_t volumes_[2]; //0 - right; 1 - left;

	enum ansys_face_type
	{
		mixed = 0,
		linear = 2,
		triangular = 3,
		quadrilateral = 4
	};

	void private_check_type_id(size_t type_id) const;
public:
	psstSimpleFace();
	psstSimpleFace(size_t type);
	psstSimpleFace(const psstSimpleFace &);
	psstSimpleFace(const size_t * vertex_index, size_t vol_left, size_t vol_right, size_t type);
	psstSimpleFace(const std::vector<size_t> & vertex_index, size_t vol_left, size_t vol_right, size_t type);

	virtual void initialize();
	virtual psstInterfaceSimpleFace * clone() const;
	virtual void copy(const psstSimpleFace &);

	virtual size_t type() const { return type_; };
	virtual void set_type(size_t type_id);

	virtual size_t number_of_vertices() const { return vertex_number_; };
	virtual size_t & vertex(size_t i) { return const_cast<size_t &>(static_cast<const psstSimpleFace &>(*this).vertex(i)); };
	virtual const size_t & vertex(size_t i) const;

	virtual size_t & cell(size_t i) { return const_cast<size_t &>(static_cast<const psstSimpleFace &>(*this).cell(i)); }
	virtual const size_t & cell(size_t i) const;

	virtual size_t & zone_id() { return zone_id_; };
	virtual const size_t & zone_id() const { return zone_id_; };

	psstSimpleFace & operator=(const psstSimpleFace & face) { copy(face); return *this; };
	virtual ~psstSimpleFace();
};

class psstSimpleCell : public psstInterfaceSimpleCell
{
private:
	size_t type_;
	size_t zone_id_;

	enum ansys_cell_type
	{
		mixed = 0,	
		triangular = 1, // 3 nodes per cell; 3 faces per cell
		tetrahedral = 2,	// 4 nodes per cell; 4 faces per cell
		quadrilateral = 3,	// 4 nodes per cell; 4 faces per cell
		hexahedral = 4, // 8 nodes per cell; 6 faces per cell
		pyramid = 5, // 5 nodes per cell; 5 faces per cell
		wedge = 6 // 6 nodes per cell; 5 faces per cell
	};
public:

	psstSimpleCell();
	psstSimpleCell(const psstSimpleCell &);
	psstSimpleCell(size_t type);
	
	virtual psstInterfaceSimpleCell * clone() const;
	virtual void copy(const psstSimpleCell &);

	virtual size_t & zone_id() { return zone_id_; };
	virtual const size_t & zone_id() const { return zone_id_; };

	virtual size_t type() const { return type_; };
	virtual void set_type(size_t type_id) { type_ = type_id; };
	virtual ~psstSimpleCell() {};
};

class psstSimpleStorage : public psstInterfaceSimpleStorage
{
private:
	size_t zone_id_;
	size_t first_index_;
	size_t last_index_;
	size_t boundary_id_;
	size_t storage_size_;
	size_t children_zone_id_;
	size_t parent_zone_id_;
	std::string zone_name_;

	enum ansys_bounary_type
	{
		interior_cell = 1,
		interior_face = 2,
		wall = 3,
		pressure_inlet = 4,
		pressure_outlet = 5,
		symmetry = 7,
		periodic_shadow = 8,
		pressure_far_field = 9,
		velocity_inlet = 10,
		periodic = 12,
		radiator = 14,
		mass_flow_inlet = 20,
		interface_ = 24,
		parent_face = 31,
		parent_cell = 32,
		outflow = 36,
		axis = 37
	};

public:
	psstSimpleStorage();
	psstSimpleStorage(const psstSimpleStorage &);
	psstSimpleStorage(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type);
	psstSimpleStorage(size_t zone_id, const std::string & name);

	virtual psstInterfaceSimpleStorage * clone() const;
	virtual void copy(const psstSimpleStorage &);

	virtual size_t zone_id() const { return zone_id_; };
	virtual void set_zone_id(size_t zone_id) { zone_id_ = zone_id; };

	virtual size_t first_index() const { return first_index_; };
	virtual size_t last_index() const { return last_index_; };
	virtual void set_range(size_t first_index, size_t last_index);
	virtual size_t storage_size() const { return storage_size_; };

	virtual size_t parent_zone_id() const { return parent_zone_id_; };
	virtual void set_parent_zone_id(size_t zone_id) { parent_zone_id_ = zone_id; };

	virtual size_t children_zone_id() const { return children_zone_id_; };
	virtual void set_children_zone_id(size_t zone_id) { children_zone_id_ = zone_id; };

	virtual size_t boundary_type() const { return boundary_id_; };
	virtual void set_boundary_type(size_t boundary_type) { boundary_id_ = boundary_type; };
	virtual bool is_interior_face() const { return boundary_id_ == interior_face; };
	virtual bool is_interior_cell() const { return boundary_id_ == interior_cell; };
	virtual bool is_parent_face() const { return boundary_id_ == parent_face; };
	virtual bool is_parent_cell() const { return boundary_id_ == parent_cell; };
	virtual bool is_boundary() const;
	virtual const std::string & zone_name() const { return zone_name_; };
	virtual void set_zone_name(const std::string & name) { zone_name_ = name; };

	virtual ~psstSimpleStorage() {};
};

class psstSimpleConnection : public psstInterfaceSimpleConnection
{
private:
	size_t parent_index_;
	size_t number_of_children_;
	size_t * children_;
public:
	psstSimpleConnection();
	psstSimpleConnection(const psstSimpleConnection &);
	psstSimpleConnection(size_t * children, size_t number_of_children, size_t parent_index);

	virtual psstInterfaceSimpleConnection * clone() const;
	virtual size_t parent_index() const { return parent_index_; };
	virtual size_t & child(size_t i) { return const_cast<size_t &>(static_cast<const psstSimpleConnection &>(*this).child(i)); };
	virtual const size_t & child(size_t i) const;
	virtual size_t number_of_children() const { return number_of_children_; };

	virtual psstSimpleConnection & operator=(const psstSimpleConnection &);
	virtual ~psstSimpleConnection();
};

class psstAnsysInputGrid : public psstInterfaceInputGrid
{
private:
	size_t dimension_;
	std::vector<psstSimpleVertex> vertices_;

	std::vector<psstSimpleFace> faces_;
	std::vector<psstSimpleConnection> faces_connections_;
	
	std::vector<psstSimpleCell> cells_;
	std::vector<psstSimpleConnection> cells_connections_;

	std::vector<psstSimpleStorage> zones_storage_;
	std::map<size_t, size_t> zones_storage_navigator_;

	bool private_is_zone_storage_exist(size_t zone_id) const;
public:
	virtual void initialize();
	virtual psstInterfaceInputGrid * clone() const;
	virtual size_t dimension() const { return dimension_; };
	virtual void set_dimension(size_t i);

	virtual void reserve_vertices(size_t vert_numb) { vertices_.reserve(vert_numb); };
	virtual size_t vertices_number() const { return vertices_.size(); };
	virtual void add_vertex(const double *, size_t);
	virtual void add_vertex(const std::vector<double> &);
	virtual void add_vertex(const psstInterfaceSimpleVertex & vertex);

	virtual void reserve_faces(size_t face_numb) { faces_.reserve(face_numb); };
	virtual size_t faces_number() const { return faces_.size(); };
	virtual void add_face(const size_t * vertex_index, size_t vol_left, size_t vol_right, size_t type, size_t zone_id);
	virtual void add_face(const std::vector<size_t> & vertex_index, size_t vol_left, size_t vol_right, size_t type, size_t zone_id);
	virtual void add_face(const psstInterfaceSimpleFace & face);
	virtual void add_face_connection(size_t * children, size_t number_of_children, size_t parent_index);

	virtual void reserve_cells(size_t cell_numb) { cells_.reserve(cell_numb); };
	virtual size_t cells_number() const { return cells_.size(); };
	virtual void add_cell(size_t type, size_t zone_id);
	virtual void add_cell(const psstInterfaceSimpleCell & cell);
	virtual void add_cell_connection(size_t * children, size_t number_of_children, size_t parent_index);

	virtual size_t number_of_zones() const { return zones_storage_.size(); };
	virtual void add_zones_storage(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type);
	virtual void add_zones_storage(size_t zone_id, const std::string & name);
	virtual void add_zones_storage(const psstInterfaceSimpleStorage & zone_storage);
	virtual void add_zones_storage_connection(size_t first_index, size_t last_index, size_t parent_zone_id, size_t children_zone_id);

	virtual psstGridFV * create_finite_volume_grid() const;

	virtual ~psstAnsysInputGrid() {};
};

class psstAnsysBuilder : public psstInterfaceGridBuilder
{
private:

	psstInterfaceInputGrid * grid_;
	
	void private_grid_ptr_check_() const;
public:
	virtual void reset();
	virtual psstInterfaceInputGrid * get_grid() const;
	virtual void read_dimension(size_t dim) const;
	virtual void read_vertices(size_t zone_id, size_t start_index, size_t end_index, size_t dimension, std::istream * input_stream) const;
	virtual void read_faces(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type, size_t face_type, std::istream * input_stream) const;
	virtual void read_cells(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type, size_t cell_type, std::istream * input_stream) const;
	virtual void read_zones(size_t zone_id, const std::string zone_name, std::istream * input_stream) const;
	virtual void read_face_tree(size_t first_index, size_t last_index, size_t parent_id, size_t children_id, std::istream * input_stream) const;
	virtual void read_cell_tree(size_t first_index, size_t last_index, size_t parent_id, size_t children_id, std::istream * input_stream) const;
	virtual ~psstAnsysBuilder() { delete grid_; };
};

class psstAnsysDirector : public psstInterfaceGridDirector
{
	psstInterfaceGridBuilder * builder_;
	std::ifstream * in_file_;

	size_t operation_id_;
	size_t dimension_;
	size_t first_index_;
	size_t last_index_;
	size_t zone_id_;
	size_t face_type_;
	size_t bc_type_;
	size_t parent_id_;
	size_t children_id_;
	std::string zone_name_;

	size_t line_pos_;
	std::string current_line_;

	size_t private_parse_number(size_t numb);
	std::string private_parse_name();

	enum ansys_operation_id
	{
		read_operation_id = 0,
		read_dimension = 2,
		read_vertices = 10,
		read_cells = 12,
		read_faces = 13,
		read_zone_39 = 39,
		read_zone_45 = 45,
		read_cell_tree = 58,
		read_face_tree = 59
	};

public:
	virtual void set_builder(psstInterfaceGridBuilder * builder);
	virtual size_t operation_id() const { return operation_id_; };
	virtual void set_operation_id(size_t operation_id) { operation_id_ = operation_id; };
	virtual void build_step();
	virtual psstInterfaceInputGrid * construct_grid(std::ifstream &);
	virtual ~psstAnsysDirector() {};
};

#endif
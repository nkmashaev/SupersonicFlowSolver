#include <vector>
#include <map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <string>
#include <iomanip>
#include "psstAnsysGridReader.h"
#include "psstFiniteVolumeGrid.h"
#include "psstFiniteVolumeGridTools.h"

size_t get_number_of_vertices(size_t type_id)
{
	switch (type_id)
	{
	case 0: //mixed
		return 0;
	case 2: //linear
		return 2;
	case 3: //triangular
		return 3;
	case 4:	//quadrilateral
		return 4;
	default:
		throw std::invalid_argument("Error: Unknown face type!");
	}

	return 0;
}

void check_dimension(size_t dim)
{
	if (dim > 3)
		throw std::length_error("Error: Invalid dimension!");
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstSimpleVertex																																													|
|																																																		|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstSimpleVertex::psstSimpleVertex(size_t dim) :
	x_(0),
	dimension_(0)
{
	private_set_dimension(dim);
	x_ = new double[dimension()];
	initialize();
}

psstSimpleVertex::psstSimpleVertex(const psstSimpleVertex & vertex) :
	x_(0),
	dimension_(0)
{
	private_set_dimension(vertex.dimension());
	x_ = new double[dimension()];
	
	for (size_t i = 0; i < dimension(); ++i)
		x_[i] = vertex.component(i);
}

psstSimpleVertex::psstSimpleVertex(const double * components, size_t dim) :
	x_(0),
	dimension_(0)
{
	if (components == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");

	private_set_dimension(dim);
	x_ = new double[dimension()];

	for (size_t i = 0; i < dimension(); ++i)
		x_[i] = components[i];
}

psstSimpleVertex::psstSimpleVertex(const std::vector<double> & components) :
	x_(0),
	dimension_(0)
{
	private_set_dimension(components.size());
	x_ = new double[dimension()];

	for (size_t i = 0; i < dimension(); ++i)
		x_[i] = components[i];
}


void psstSimpleVertex::initialize()
{
	for (size_t i = 0; i < dimension(); ++i)
		component(i) = 0.0;
}

psstInterfaceSimpleVertex * psstSimpleVertex::clone() const
{
	psstInterfaceSimpleVertex * new_simple_vertex = new psstSimpleVertex(*this);
	return new_simple_vertex;
}

void psstSimpleVertex::copy(const psstSimpleVertex & simple_vertex)
{
	if (&simple_vertex == this)
		return;
	delete[] x_;
	
	private_set_dimension(simple_vertex.dimension());
	x_ = new double[dimension()];

	for (size_t i = 0; i < dimension(); ++i)
		x_[i] = simple_vertex.component(i);
}

const double & psstSimpleVertex::component(size_t i) const
{
	if (i >= dimension())
		throw std::out_of_range("Error: Out of range!");
	return x_[i];
}

psstSimpleVertex & psstSimpleVertex::operator=(const psstSimpleVertex & vert)
{
	copy(vert);
	return *this;
}

void psstSimpleVertex::private_set_dimension(size_t dim)
{
	check_dimension(dim);
	dimension_ = dim;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstSimpleFace																																														|
|																																																		|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstSimpleFace::psstSimpleFace() :
	zone_id_(0),
	type_(mixed),
	vertex_number_(0),
	vertex_index_(0),
	volumes_{0, 0}
{}

psstSimpleFace::psstSimpleFace(size_t type) :
	zone_id_(0),
	type_(mixed),
	vertex_number_(0),
	vertex_index_(0),
	volumes_{ 0, 0 }
{
	set_type(type);
}

psstSimpleFace::psstSimpleFace(const psstSimpleFace & face) :
	zone_id_(0),
	type_(mixed),
	vertex_number_(0),
	vertex_index_(0),
	volumes_{ 0, 0 }
{
	vertex_index_ = new size_t[face.number_of_vertices()];

	vertex_number_ = face.number_of_vertices();
	type_ = face.type();
	zone_id() = face.zone_id();
	for (size_t i = 0; i < number_of_vertices(); ++i)
		vertex_index_[i] = face.vertex(i);
	for (size_t i = 0; i < 2; ++i)
		volumes_[i] = face.cell(i);
}

psstSimpleFace::psstSimpleFace(const size_t * vertex_index, size_t vol_left, size_t vol_right, size_t type) :
	zone_id_(0),
	type_(mixed),
	vertex_number_(0),
	vertex_index_(0),
	volumes_{ 0, 0 }
{
	if (vertex_index == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");

	set_type(type);
	for (size_t i = 0; i < number_of_vertices(); ++i)
	{
		vertex(i) = vertex_index[i];
	}
	cell(0) = vol_right;
	cell(1) = vol_left;
}

psstSimpleFace::psstSimpleFace(const std::vector<size_t> & vertex_index, size_t vol_left, size_t vol_right, size_t type) :
	zone_id_(0),
	type_(mixed),
	vertex_number_(0),
	vertex_index_(0),
	volumes_{ 0, 0 }
{
	set_type(type);
	if (number_of_vertices() != vertex_index.size())
		throw std::invalid_argument("Error: Invalid number of vertex!");

	for (size_t i = 0; i < number_of_vertices(); ++i)
		vertex(i) = vertex_index[i];
	cell(0) = vol_right;
	cell(1) = vol_left;
}

void psstSimpleFace::initialize()
{
	for (size_t i = 0; i < number_of_vertices(); ++i)
	{
		vertex(i) = 0;
	}
	cell(0) = 0;
	cell(1) = 0;
}

psstInterfaceSimpleFace * psstSimpleFace::clone() const
{
	psstInterfaceSimpleFace * new_simple_face = new psstSimpleFace(*this);
	return new_simple_face;
}

void psstSimpleFace::copy(const psstSimpleFace & face)
{
	if (&face == this)
		return;

	delete[] vertex_index_;
	set_type(face.type());
	zone_id() = face.zone_id();
	for (size_t i = 0; i < number_of_vertices(); ++i)
	{
		vertex(i) = face.vertex(i);
	}

	for (size_t i = 0; i < 2; ++i)
	{
		cell(i) = face.cell(i);
	}
}

void psstSimpleFace::set_type(size_t type_id)
{
	switch (type_id)
	{
	case mixed:
		type_ = mixed;
		vertex_number_ = 0;
		break;
	case linear:
		type_ = linear;
		vertex_number_ = 2;
		break;
	case triangular:
		type_ = triangular;
		vertex_number_ = 3;
		break;
	case quadrilateral:
		type_ = quadrilateral;
		vertex_number_ = 4;
		break;
	default:
		throw std::invalid_argument("Error: Unknown face type!");
	}

	vertex_index_ = new size_t[number_of_vertices()];
	initialize();
}

const size_t & psstSimpleFace::vertex(size_t i) const
{
	if (i > number_of_vertices())
		throw std::out_of_range("Error: Out of range!");
	return vertex_index_[i];
}

const size_t & psstSimpleFace::cell(size_t i) const
{
	if (i > 1)
		throw std::out_of_range("Error: Out of range!");
	return volumes_[i];
}

psstSimpleFace::~psstSimpleFace()
{
	delete[]vertex_index_;
}

void psstSimpleFace::private_check_type_id(size_t type_id) const
{
	switch (type_id)
	{
	case mixed :
		return;
	case linear:
		return;
	case triangular :
		return;
	case quadrilateral :
		return;
	default:
		throw std::invalid_argument("Error: Unknown face type!");
	}
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstSimpleCell																																														|
|																																																		|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstSimpleCell::psstSimpleCell() :
	type_(0),
	zone_id_(0)
{}

psstSimpleCell::psstSimpleCell(const psstSimpleCell & cell) :
	type_(cell.type_),
	zone_id_(cell.zone_id_)
{}

psstSimpleCell::psstSimpleCell(size_t type) :
	type_(type),
	zone_id_(0)
{}

psstInterfaceSimpleCell * psstSimpleCell::clone() const
{
	psstInterfaceSimpleCell * new_cell = new psstSimpleCell(*this);
	return new_cell;
}

void psstSimpleCell::copy(const psstSimpleCell & cell)
{
	if (this == &cell)
		return;
	set_type(cell.type());
	zone_id() = cell.zone_id();
}


/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstSimpleStorage																																													|
|																																																		|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
psstSimpleStorage::psstSimpleStorage() :
	zone_id_(0),
	first_index_(0),
	last_index_(0),
	boundary_id_(0),
	storage_size_(0),
	children_zone_id_(0),
	parent_zone_id_(0),
	zone_name_()
{}

psstSimpleStorage::psstSimpleStorage(const psstSimpleStorage & storage) :
	zone_id_(storage.zone_id_),
	first_index_(storage.first_index_),
	last_index_(storage.last_index_),
	boundary_id_(storage.boundary_id_),
	storage_size_(storage.storage_size_),
	children_zone_id_(storage.children_zone_id_),
	parent_zone_id_(storage.parent_zone_id_),
	zone_name_(storage.zone_name_)
{}

psstSimpleStorage::psstSimpleStorage(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type) :
	zone_id_(0),
	first_index_(0),
	last_index_(0),
	boundary_id_(0),
	storage_size_(0),
	children_zone_id_(0),
	parent_zone_id_(0),
	zone_name_()
{
	set_range(first_index, last_index);
	zone_id_ = zone_id;
	boundary_id_ = bc_type;
}

psstSimpleStorage::psstSimpleStorage(size_t zone_id, const std::string & name) :
	zone_id_(zone_id),
	first_index_(0),
	last_index_(0),
	boundary_id_(0),
	storage_size_(0),
	children_zone_id_(0),
	parent_zone_id_(0),
	zone_name_(name)
{}

psstInterfaceSimpleStorage * psstSimpleStorage::clone() const
{
	psstInterfaceSimpleStorage * new_storage = new psstSimpleStorage(*this);
	return new_storage;
}

void psstSimpleStorage::set_range(size_t first_index, size_t last_index)
{
	if (last_index < first_index)
		throw std::invalid_argument("Error: Expected first index to be less than the last!");
	first_index_ = first_index;
	last_index_ = last_index;
	storage_size_ = last_index_ - first_index_ + 1;
}

void psstSimpleStorage::copy(const psstSimpleStorage & storage)
{
	if (&storage == this)
		return;

	set_range(storage.first_index(), storage.last_index());
	set_zone_name(storage.zone_name());
	set_boundary_type(storage.boundary_type());
	set_zone_id(storage.zone_id());
}

bool psstSimpleStorage::is_boundary() const
{
	return (!is_interior_cell() && !is_interior_face());
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstSimpleConnection																																												|
|																																																		|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstSimpleConnection::psstSimpleConnection() :
	parent_index_(0),
	number_of_children_(0),
	children_(0)
{}

psstSimpleConnection::psstSimpleConnection(const psstSimpleConnection & connection) :
	parent_index_(0),
	number_of_children_(0),
	children_(0)
{
	number_of_children_ = connection.number_of_children();
	parent_index_ = connection.parent_index();
	
	children_ = new size_t[number_of_children()];
	for (size_t i = 0; i < number_of_children(); ++i)
		child(i) = connection.child(i);
}

psstSimpleConnection::psstSimpleConnection(size_t * children, size_t number_of_children, size_t parent_index)
{
	if (children == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");

	number_of_children_ = number_of_children;
	parent_index_ = parent_index;
	children_ = new size_t[psstSimpleConnection::number_of_children()];
	for (size_t i = 0; i < psstSimpleConnection::number_of_children(); ++i)
		child(i) = children[i];
}

psstInterfaceSimpleConnection * psstSimpleConnection::clone() const
{
	psstInterfaceSimpleConnection * new_connection = new psstSimpleConnection(*this);
	return new_connection;
}

const size_t & psstSimpleConnection::child(size_t i) const
{
	if (i > number_of_children())
		throw std::out_of_range("Error: Out of range!");
	return children_[i];
}

psstSimpleConnection & psstSimpleConnection::operator=(const psstSimpleConnection & connection)
{
	if (&connection == this)
		return *this;

	delete[] children_;
	number_of_children_ = connection.number_of_children();
	parent_index_ = connection.parent_index();
	for (size_t i = 0; i < number_of_children(); ++i)
		child(i) = connection.child(i);
	return *this;
}

psstSimpleConnection::~psstSimpleConnection()
{
	delete[] children_;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceInputGrid																																												|
|	psstAnsysInputGrid																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstAnsysInputGrid::initialize()
{

}

psstInterfaceInputGrid * psstAnsysInputGrid::clone() const
{
	psstInterfaceInputGrid * new_grid = new psstAnsysInputGrid(*this);
	return new_grid;
}

void psstAnsysInputGrid::set_dimension(size_t dim)
{
	check_dimension(dim);
	dimension_ = dim;
}

void psstAnsysInputGrid::add_vertex(const double * components, size_t dim)
{
	if (components == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");
	psstSimpleVertex simple_vertex(components, dim);
	vertices_.push_back(simple_vertex);
}

void psstAnsysInputGrid::add_vertex(const std::vector<double> & components)
{
	psstSimpleVertex simple_vertex(components);
	vertices_.push_back(simple_vertex);
}

void psstAnsysInputGrid::add_vertex(const psstInterfaceSimpleVertex & vertex)
{
	const psstSimpleVertex * vertex1 = dynamic_cast<const psstSimpleVertex *>(&vertex);
	if (vertex1 == 0)
		throw std::bad_cast();
	vertices_.push_back(*vertex1);
}

void psstAnsysInputGrid::add_face(const size_t * vertex_index, size_t vol_left, size_t vol_right, size_t type, size_t zone_id)
{
	if (vertex_index == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");
	psstSimpleFace simple_face(vertex_index, vol_left, vol_right, type);
	simple_face.zone_id() = zone_id;
	faces_.push_back(simple_face);
}

void psstAnsysInputGrid::add_face(const std::vector<size_t> & vertex_index, size_t vol_left, size_t vol_right, size_t type, size_t zone_id)
{
	psstSimpleFace simple_face(vertex_index, vol_left, vol_right, type);
	simple_face.zone_id() = zone_id;
	faces_.push_back(simple_face);
}

void psstAnsysInputGrid::add_face(const psstInterfaceSimpleFace & face)
{
	const psstSimpleFace * face1 = dynamic_cast<const psstSimpleFace *>(&face);
	if (face1 == 0)
		throw std::bad_cast();
	faces_.push_back(*face1);
}

void psstAnsysInputGrid::add_face_connection(size_t * children, size_t number_of_children, size_t parent_index)
{
	psstSimpleConnection simple_connection(children, number_of_children, parent_index);
	faces_connections_.push_back(simple_connection);
}

void psstAnsysInputGrid::add_cell(size_t type, size_t zone_id)
{
	psstSimpleCell cell(type);
	cell.zone_id() = zone_id;
	cells_.push_back(cell);
}

void psstAnsysInputGrid::add_cell(const psstInterfaceSimpleCell & cell)
{
	const psstSimpleCell * cell1 = dynamic_cast<const psstSimpleCell *>(&cell);
	if (cell1 == 0)
		throw std::bad_cast();
	cells_.push_back(*cell1);
}

void psstAnsysInputGrid::add_cell_connection(size_t * children, size_t number_of_children, size_t parent_index)
{
	psstSimpleConnection simple_connection(children, number_of_children, parent_index);
	cells_connections_.push_back(simple_connection);
}

void psstAnsysInputGrid::add_zones_storage(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type)
{
	if (!private_is_zone_storage_exist(zone_id))
	{
		size_t new_storage_index = zones_storage_.size();
		psstSimpleStorage new_storage(zone_id, first_index, last_index, bc_type);
		zones_storage_.push_back(new_storage);
		zones_storage_navigator_[zone_id] = new_storage_index;
	}
	else
	{
		size_t index = zones_storage_navigator_.find(zone_id)->second;
		zones_storage_[index].set_range(first_index, last_index);
		zones_storage_[index].set_zone_id(zone_id);
		zones_storage_[index].set_boundary_type(bc_type);
	}
}

void psstAnsysInputGrid::add_zones_storage(size_t zone_id, const std::string & name)
{
	if (!private_is_zone_storage_exist(zone_id))
	{
		size_t new_storage_index = zones_storage_.size();
		psstSimpleStorage new_storage(zone_id, name);
		zones_storage_.push_back(new_storage);
		zones_storage_navigator_[zone_id] = new_storage_index;
	}
	else
	{
		size_t index = zones_storage_navigator_.find(zone_id)->second;
		zones_storage_[index].set_zone_id(zone_id);
		zones_storage_[index].set_zone_name(name);
	}
}

void psstAnsysInputGrid::add_zones_storage(const psstInterfaceSimpleStorage & zone_storage)
{
	const psstSimpleStorage * zone_storage1 = dynamic_cast<const psstSimpleStorage *>(&zone_storage);
	if (zone_storage1 == 0)
		throw std::bad_cast();

	size_t zone_id = zone_storage1->zone_id();
	if (!private_is_zone_storage_exist(zone_id))
	{
		size_t new_storage_index = zones_storage_.size();
		zones_storage_.push_back(*zone_storage1);
		zones_storage_navigator_[zone_id] = new_storage_index;
	}
	else
	{
		size_t index = zones_storage_navigator_.find(zone_id)->second;
		zones_storage_[index] = *zone_storage1;
	}
}

void psstAnsysInputGrid::add_zones_storage_connection(size_t first_index, size_t last_index, size_t parent_zone_id, size_t children_zone_id)
{
	size_t zone_id = parent_zone_id;
	if (!private_is_zone_storage_exist(zone_id))
	{
		size_t new_storage_index = zones_storage_.size();
		psstSimpleStorage new_storage(parent_zone_id, first_index, last_index, 31);
		new_storage.set_children_zone_id(children_zone_id);
		zones_storage_.push_back(new_storage);
		zones_storage_navigator_[zone_id] = new_storage_index;
	}
	else
	{
		size_t index = zones_storage_navigator_.find(zone_id)->second;
		zones_storage_[index].set_children_zone_id(children_zone_id);
	}

	zone_id = children_zone_id;
	if (!private_is_zone_storage_exist(zone_id))
	{
		size_t new_storage_index = zones_storage_.size();
		psstSimpleStorage new_storage;
		new_storage.set_parent_zone_id(parent_zone_id);
		zones_storage_.push_back(new_storage);
		zones_storage_navigator_[zone_id] = new_storage_index;
	}
	else
	{
		size_t index = zones_storage_navigator_.find(zone_id)->second;
		zones_storage_[index].set_parent_zone_id(parent_zone_id);
	}
}

psstGridFV * psstAnsysInputGrid::create_finite_volume_grid() const
{
	psstGridFV * finite_volume_grid = new psstGridFV;

	//Задание размерности
	size_t dim = dimension();
	finite_volume_grid->set_dimension(dim);

	//Считывание координат
	std::vector<std::shared_ptr<psstInterfaceVertex>> * vs_ptr = finite_volume_grid->vertex_storage();
	vs_ptr->resize(vertices_number());
	for (size_t i = 0; i < vs_ptr->size(); ++i)
	{
		vs_ptr->at(i) = std::make_shared<psstVertex>(dim);
		for (size_t j = 0; j < dim; ++j)
		{
			vs_ptr->at(i)->component(j) = vertices_[i].component(j);
		}
		vs_ptr->at(i)->id() = i + 1;
	}

	//Выделение памяти под грани и объемы
	std::map<size_t, psstFaceStorage> * fs_ptr = finite_volume_grid->face_storage();
	std::map<size_t, psstVolumeStorage> * cs_ptr = finite_volume_grid->volume_storage();
	for (size_t i = 0; i < number_of_zones(); ++i)
	{
		size_t zone_id = zones_storage_[i].zone_id();
		size_t bc_type = zones_storage_[i].boundary_type();
		if (bc_type == 1 || bc_type == 32) //interior fluent
		{
			cs_ptr->operator[](zone_id).volumes_.resize(zones_storage_[i].storage_size());
			cs_ptr->operator[](zone_id).parent_zone_id_ = zones_storage_[i].parent_zone_id();
			cs_ptr->operator[](zone_id).children_zone_id_ = zones_storage_[i].children_zone_id();
			cs_ptr->operator[](zone_id).bc_type_id_ = zones_storage_[i].boundary_type();
			cs_ptr->operator[](zone_id).name_ = zones_storage_[i].zone_name();
			cs_ptr->operator[](zone_id).first_index_ = zones_storage_[i].first_index();
			cs_ptr->operator[](zone_id).last_index_ = zones_storage_[i].last_index();

			for (size_t j = 0; j < cs_ptr->at(zone_id).volumes_.size(); ++j)
			{
				size_t loc = zones_storage_[i].first_index() + j - 1;
				cs_ptr->at(zone_id).volumes_.at(j) = std::make_shared<psstVolume>(dim);
				cs_ptr->at(zone_id).volumes_.at(j)->set_type(cells_[loc].type());
				cs_ptr->at(zone_id).volumes_.at(j)->id() = loc + 1;
			}
		}
		else if (bc_type == 2 || bc_type == 31)
		{
			fs_ptr->operator[](zone_id).faces_.resize(zones_storage_[i].storage_size());
			fs_ptr->operator[](zone_id).parent_zone_id_ = zones_storage_[i].parent_zone_id();
			fs_ptr->operator[](zone_id).children_zone_id_ = zones_storage_[i].children_zone_id();
			fs_ptr->operator[](zone_id).bc_type_id_ = zones_storage_[i].boundary_type();
			fs_ptr->operator[](zone_id).name_ = zones_storage_[i].zone_name();
			fs_ptr->operator[](zone_id).first_index_ = zones_storage_[i].first_index();
			fs_ptr->operator[](zone_id).last_index_ = zones_storage_[i].last_index();

			for (size_t j = 0; j < fs_ptr->at(zone_id).faces_.size(); ++j)
			{
				size_t loc = zones_storage_[i].first_index() + j - 1;
				fs_ptr->at(zone_id).faces_.at(j) = std::make_shared<psstFace>(dim);
				fs_ptr->at(zone_id).faces_.at(j)->set_type(faces_[loc].type());
				fs_ptr->at(zone_id).faces_.at(j)->id() = loc + 1;
			}
		}
		else
		{
			cs_ptr->operator[](zone_id).volumes_.resize(zones_storage_[i].storage_size());
			cs_ptr->operator[](zone_id).parent_zone_id_ = zones_storage_[i].parent_zone_id();
			cs_ptr->operator[](zone_id).children_zone_id_ = zones_storage_[i].children_zone_id();
			cs_ptr->operator[](zone_id).bc_type_id_ = zones_storage_[i].boundary_type();
			cs_ptr->operator[](zone_id).name_ = zones_storage_[i].zone_name();
			cs_ptr->operator[](zone_id).first_index_ = zones_storage_[i].first_index();
			cs_ptr->operator[](zone_id).last_index_ = zones_storage_[i].last_index();

			fs_ptr->operator[](zone_id).faces_.resize(zones_storage_[i].storage_size());
			fs_ptr->operator[](zone_id).parent_zone_id_ = zones_storage_[i].parent_zone_id();
			fs_ptr->operator[](zone_id).children_zone_id_ = zones_storage_[i].children_zone_id();
			fs_ptr->operator[](zone_id).bc_type_id_ = zones_storage_[i].boundary_type();
			fs_ptr->operator[](zone_id).name_ = zones_storage_[i].zone_name();
			fs_ptr->operator[](zone_id).first_index_ = zones_storage_[i].first_index();
			fs_ptr->operator[](zone_id).last_index_ = zones_storage_[i].last_index();

			for (size_t j = 0; j < fs_ptr->at(zone_id).faces_.size(); ++j)
			{
				size_t loc = zones_storage_[i].first_index() + j - 1;
				fs_ptr->at(zone_id).faces_.at(j) = std::make_shared<psstFace>(dim);
				fs_ptr->at(zone_id).faces_.at(j)->set_boundary_status(true);
				fs_ptr->at(zone_id).faces_.at(j)->set_type(faces_[loc].type());
				fs_ptr->at(zone_id).faces_.at(j)->id() = loc + 1;

				cs_ptr->at(zone_id).volumes_.at(j) = std::make_shared<psstVolume>(dim);
				cs_ptr->at(zone_id).volumes_.at(j)->set_boundary_status(true);
				cs_ptr->at(zone_id).volumes_.at(j)->id() = 0;
			}
		}
	}
	finite_volume_grid->set_number_of_faces(faces_.size());
	finite_volume_grid->set_number_of_inner_volumes(cells_.size());

	//Привязка отношений детей-ребенок
	for (size_t i = 0; i < faces_connections_.size(); ++i)
	{
		size_t parent_zone_id = faces_[faces_connections_[i].parent_index() - 1].zone_id();
		size_t zone_id = zones_storage_navigator_.find(parent_zone_id)->second;
		size_t first_index = zones_storage_.at(zone_id).first_index();
		size_t parent_local_index = faces_connections_[i].parent_index() - first_index;
			
		for (size_t j = 0; j < faces_connections_[i].number_of_children(); ++j)
		{
			size_t child_zone_id = faces_[faces_connections_[i].child(j) - 1].zone_id();
			size_t curr_zone_index = zones_storage_navigator_.find(child_zone_id)->second;
			size_t child_local_index = faces_connections_[i].child(j) - zones_storage_.at(curr_zone_index).first_index();

			fs_ptr->at(parent_zone_id).faces_[parent_local_index]->add_child(fs_ptr->at(child_zone_id).faces_[child_local_index]);
			fs_ptr->at(child_zone_id).faces_[child_local_index]->parent() = fs_ptr->at(parent_zone_id).faces_[parent_local_index];
		}
	}

	for (size_t i = 0; i < cells_connections_.size(); ++i)
	{
		size_t parent_zone_id = cells_[cells_connections_[i].parent_index() - 1].zone_id();
		size_t zone_id = zones_storage_navigator_.find(parent_zone_id)->second;
		size_t first_index = zones_storage_.at(zone_id).first_index();
		size_t parent_local_index = cells_connections_[i].parent_index() - first_index;

		for (size_t j = 0; j < cells_connections_[i].number_of_children(); ++j)
		{
			size_t child_zone_id = cells_[cells_connections_[i].child(j) - 1].zone_id();
			size_t curr_zone_index = zones_storage_navigator_.find(child_zone_id)->second;
			size_t child_local_index = cells_connections_[i].child(j) - zones_storage_.at(curr_zone_index).first_index();

			cs_ptr->at(parent_zone_id).volumes_[parent_local_index]->add_child (cs_ptr->at(child_zone_id).volumes_[child_local_index]);
			cs_ptr->at(child_zone_id).volumes_[child_local_index]->parent() = cs_ptr->at(parent_zone_id).volumes_[parent_local_index];
		}
	}

	//Привязка граней к объемам
	for (size_t i = 0; i < number_of_zones(); ++i)
	{
		size_t zone_id = zones_storage_[i].zone_id();
		size_t bc_type = zones_storage_[i].boundary_type();
		if (bc_type != 1 && bc_type != 32)
		{
			size_t loc = 0;
			for (size_t j = 0; j < fs_ptr->at(zone_id).faces_.size(); ++j)
			{
				loc = zones_storage_[i].first_index() + j - 1;

				for (size_t k = 0; k < faces_[loc].number_of_vertices(); ++k)
				{
					fs_ptr->at(zone_id).faces_.at(j)->add_vertex(vs_ptr->at(faces_[loc].vertex(k) - 1));
				}

				for (size_t k = 0; k < 2; ++k)
				{
					size_t neighbour_cell = faces_.at(loc).cell(k);
					size_t vol_zone_id = 0;
					if (bc_type != 31)
					{
						if (neighbour_cell == 0)
						{
							neighbour_cell = j;
							vol_zone_id = zone_id;
						}
						else
						{
							neighbour_cell = faces_.at(loc).cell(k) - 1;
							vol_zone_id = cells_.at(neighbour_cell).zone_id();
							size_t zone_id1 = zones_storage_navigator_.find(vol_zone_id)->second;
							size_t first_index = zones_storage_.at(zone_id1).first_index();
							neighbour_cell = faces_.at(loc).cell(k) - first_index;
						}
						fs_ptr->at(zone_id).faces_.at(j)->neighbour(k) = cs_ptr->at(vol_zone_id).volumes_.at(neighbour_cell);
						cs_ptr->at(vol_zone_id).volumes_.at(neighbour_cell)->add_face(fs_ptr->at(zone_id).faces_.at(j));
					}
					else
					{
						if (neighbour_cell == 0)
						{
							fs_ptr->at(zone_id).faces_.at(j)->neighbour(k) = 0;
						}
						else
						{
							neighbour_cell = faces_.at(loc).cell(k) - 1;
							vol_zone_id = cells_.at(neighbour_cell).zone_id();
							size_t zone_id1 = zones_storage_navigator_.find(vol_zone_id)->second;
							size_t first_index = zones_storage_.at(zone_id1).first_index();
							neighbour_cell = faces_.at(loc).cell(k) - first_index;
							fs_ptr->at(zone_id).faces_.at(j)->neighbour(k) = cs_ptr->at(vol_zone_id).volumes_.at(neighbour_cell);
						}
					}
				}
				fs_ptr->at(zone_id).faces_.at(j)->initialize();
			}
		}		
	}
	
	//Инициализация
	double min = -1.0;
	for (size_t i = 0; i < number_of_zones(); ++i)
	{
		size_t zone_id = zones_storage_[i].zone_id();
		size_t bc_type = zones_storage_[i].boundary_type();
		if (bc_type != 2 && bc_type != 31 && bc_type != 32) //interior fluent
		{
			for (size_t j = 0; j < cs_ptr->at(zone_id).volumes_.size(); ++j)
			{
				cs_ptr->at(zone_id).volumes_.at(j)->initialize();
				
				if (bc_type == 1)
				{
					if (min > cs_ptr->at(zone_id).volumes_.at(j)->volume() || min < 0.0)
						min = cs_ptr->at(zone_id).volumes_.at(j)->volume();
				}
			}
		}
	}
	finite_volume_grid->set_min_vol(min);
	return finite_volume_grid;
}

bool psstAnsysInputGrid::private_is_zone_storage_exist(size_t zone_id) const
{
	std::map<size_t, size_t>::const_iterator iter = zones_storage_navigator_.find(zone_id);
	if (iter != zones_storage_navigator_.cend())
		return true;
	return false;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceGridBuilder																																											|
|	psstAnsysBuilder																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstAnsysBuilder::reset()
{
	delete grid_;
	grid_ = new psstAnsysInputGrid;
}

psstInterfaceInputGrid * psstAnsysBuilder::get_grid() const
{
	private_grid_ptr_check_();
	psstInterfaceInputGrid * grid = grid_->clone();
	return grid;
}

void psstAnsysBuilder::read_dimension(size_t dim) const
{
	std::cout << "Reading dimension section..." << std::endl;
	private_grid_ptr_check_();

	size_t dimension = dim;
	std::cout << "Dimension: " << dim << "D" << std::endl;;
	grid_->set_dimension(dim);
	std::cout << std::endl;
}

void psstAnsysBuilder::read_vertices(size_t zone_id, size_t first_index, size_t last_index, size_t dimension, std::istream * input_stream) const
{
	std::cout << "Reading vertices section..." << std::endl;
	private_grid_ptr_check_();

	size_t storage_size = last_index - first_index + 1;
	if (zone_id == 0)
	{
		grid_->reserve_vertices(storage_size);
		std::cout << "Total number of vertices: " << storage_size << std::endl;
		std::cout << "Dimension: " << dimension << "D" << std::endl;
		return;
	}

	std::cout << "zone_id = " << zone_id << ", first_index = " << first_index << ", last_index = " << last_index << std::endl;
	std::cout << "number of vertices of current zone: " << storage_size << std::endl;

	size_t dim = grid_->dimension();
	check_dimension(dim);
	double * x = new double[dim];
	for (size_t i = 0; i < storage_size; ++i)
	{
		for (size_t j = 0; j < dim; ++j)
		{
			*input_stream >> x[j];
		}
		grid_->add_vertex(x, dim);
	}
	delete[] x;
	std::cout << std::endl;
}

void psstAnsysBuilder::read_faces(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type, size_t face_type, std::istream * input_stream) const
{
	std::cout << "Reading faces section..." << std::endl;
	private_grid_ptr_check_();

	size_t storage_size = last_index - first_index + 1;
	if (zone_id == 0)
	{
		grid_->reserve_faces(storage_size);
		std::cout << "Total number of faces: " << storage_size << std::endl;
		return;
	}

	std::cout << "zone_id = " << zone_id << ", first_index = " << first_index << ", last_index = " << last_index << ", boundary_id = " << bc_type << ", face_type = " << face_type <<std::endl;
	std::cout << "number of faces of current zone: " << storage_size << std::endl;
	grid_->add_zones_storage(zone_id,first_index,last_index, bc_type);

	size_t type = face_type;
	for (size_t i = 0; i < storage_size; ++i)
	{
		if (face_type == 0)
		{
			*input_stream >> std::hex >> type;
			/*std::cout << type << " ";*/
		}
		
		
		size_t number_of_vertices = get_number_of_vertices(type);
		//std::cout << " ("<< number_of_vertices << ") ";

		//Ввод индексов вершин
		size_t * vertices = new size_t[number_of_vertices];
		for (size_t j = 0; j < number_of_vertices; ++j)
		{
			*input_stream >> std::hex >> vertices[j];
			//std::cout << vertices[j] << " ";
		}

		//Ввод индексов соседних объемов
		size_t cells[2] = {0, 0};
		for (size_t j = 0; j < 2; ++j)
		{
			*input_stream >> std::hex >> cells[j];
			//std::cout << cells[j] << " ";
		}
		//std::cout << std::endl;
		grid_->add_face(vertices, cells[1], cells[0], type, zone_id);
		delete[] vertices;
	}
	std::cout << std::endl;
}

void psstAnsysBuilder::read_cells(size_t zone_id, size_t first_index, size_t last_index, size_t bc_type, size_t cell_type, std::istream * input_stream) const
{
	std::cout << "Reading cells section..." << std::endl;
	private_grid_ptr_check_();
	
	size_t storage_size = last_index - first_index + 1;
	if (zone_id == 0)
	{
		grid_->reserve_cells(storage_size);
		std::cout << "Total number of cells: " << storage_size << std::endl;
		return;
	}

	std::cout << "zone_id = " << zone_id << ", first_index = " << first_index << ", last_index = " << last_index << ", boundary_id = " << bc_type << ", cell_type = " << cell_type << std::endl;
	std::cout << "number of cells of current zone: " << storage_size << std::endl;
	grid_->add_zones_storage(zone_id, first_index, last_index, bc_type);
	
	size_t cell_type_in = cell_type;
	for (size_t i = 0; i < storage_size; ++i)
	{
		if (cell_type == 0)
			*input_stream >> cell_type_in;
		grid_->add_cell(cell_type_in, zone_id);
	}
	std::cout << std::endl;
}

void psstAnsysBuilder::read_zones(size_t zone_id, const std::string zone_name, std::istream * input_stream) const
{
	std::cout << "Reading zones section..." << std::endl;
	private_grid_ptr_check_();
	std::cout << "zone_id = " << zone_id << ", zone_name:  " << zone_name << std::endl;
	grid_->add_zones_storage(zone_id, zone_name);
	std::cout << std::endl;
}

void psstAnsysBuilder::read_face_tree(size_t first_index, size_t last_index, size_t parent_id, size_t children_id, std::istream * input_stream) const
{
	std::cout << "Reading faces-tree section..." << std::endl;
	private_grid_ptr_check_();
	grid_->add_zones_storage_connection(first_index, last_index, parent_id, children_id);
	size_t storage_size = last_index - first_index + 1;
	
	std::cout << "first_index = " << first_index << ", last_index = " << last_index << ", parent_id = " << parent_id << ", children_id = " << children_id << std::endl;
	std::cout << "number of parent faces: " << storage_size << std::endl;

	size_t parent_index = first_index;
	size_t number_of_children = 0;
	for (size_t i = 0; i < storage_size; ++i)
	{
		parent_index = first_index + i;
		*input_stream >> std::hex >> number_of_children;
		size_t * children = new size_t[number_of_children];
		for (size_t j = 0; j < number_of_children; ++j)
			*input_stream >> std::hex >> children[j];
		grid_->add_face_connection(children, number_of_children, parent_index);
		delete[] children;
	}
	std::cout << std::endl;
}

void psstAnsysBuilder::read_cell_tree(size_t first_index, size_t last_index, size_t parent_id, size_t children_id, std::istream * input_stream) const
{
	std::cout << "Reading cells-tree section..." << std::endl;
	private_grid_ptr_check_();
	grid_->add_zones_storage_connection(first_index, last_index, parent_id, children_id);
	size_t storage_size = last_index - first_index + 1;

	std::cout << "first_index = " << first_index << ", last_index = " << last_index << ", parent_id = " << parent_id << ", children_id = " << children_id << std::endl;
	std::cout << "number of parent cells: " << storage_size << std::endl;

	size_t parent_index = first_index;
	size_t number_of_children = 0;
	for (size_t i = 0; i < storage_size; ++i)
	{
		parent_index = first_index + i;
		*input_stream >> std::hex >> number_of_children;
		size_t * children = new size_t[number_of_children];
		for (size_t j = 0; j < number_of_children; ++j)
			*input_stream >> std::hex >> children[j];
		grid_->add_cell_connection(children, number_of_children, parent_index);
		delete[] children;
	}
	std::cout << std::endl;
}

void psstAnsysBuilder::private_grid_ptr_check_() const
{
	if (grid_ == 0)
		throw std::runtime_error("Error: Current grid is not initialized!");
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceGridDirector																																											|
|	psstAnsysDirector																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstAnsysDirector::set_builder(psstInterfaceGridBuilder * builder)
{
	if (builder == 0)
		std::invalid_argument("Error: Expected builder but null_ptr found!");
	builder_ = builder;
}

void psstAnsysDirector::build_step()
{
	size_t operation_id = psstAnsysDirector::operation_id();
	size_t type = 0;
	switch (operation_id)
	{
	case read_operation_id:
		getline(*in_file_, current_line_);
		//std::cout << current_line_ << std::endl;
		line_pos_ = 0;
		set_operation_id(private_parse_number(10));
		break;
	case read_dimension:
		dimension_ = private_parse_number(16);
		builder_->read_dimension(dimension_);
		set_operation_id(read_operation_id);
		break;
	case read_vertices:

		zone_id_ = private_parse_number(16);
		first_index_ = private_parse_number(16);
		last_index_ = private_parse_number(16);
		type = private_parse_number(16);
		dimension_ = private_parse_number(16);
		if (zone_id_ != 0)
		{
			for (int i = current_line_.length(); i > 0; --i)
			{
				char ch = current_line_[i];
				if (ch == ' ' || ch == '\0')
					continue;
				else if (ch == '(')
					break;
				else
				{
					getline(*in_file_, current_line_);
					break;
				}
			}
		}
		builder_->read_vertices(zone_id_,first_index_,last_index_,dimension_, in_file_);
		set_operation_id(read_operation_id);
		break;

	case read_cells:

		zone_id_ = private_parse_number(16);
		first_index_ = private_parse_number(16);
		last_index_ = private_parse_number(16);
		bc_type_ = private_parse_number(16);
		type = private_parse_number(16);
		if (zone_id_ != 0)
		{
			for (int i = current_line_.length(); i > 0; --i)
			{
				char ch = current_line_[i];
				if (ch == ' ' || ch == '\0')
					continue;
				else if (ch == '(')
					break;
				else
				{
					getline(*in_file_, current_line_);
					break;
				}
			}
		}
		builder_->read_cells(zone_id_,first_index_,last_index_,bc_type_,type,in_file_);
		set_operation_id(read_operation_id);
		break;

	case read_faces:

		zone_id_ = private_parse_number(16);
		first_index_ = private_parse_number(16);
		last_index_ = private_parse_number(16);
		if (zone_id_ != 0)
		{
			bc_type_ = private_parse_number(16);
			face_type_ = private_parse_number(16);
			for (int i = current_line_.length(); i > 0; --i)
			{
				char ch = current_line_[i];
				if (ch == ' ' || ch == '\0')
					continue;
				else if (ch == '(')
					break;
				else
				{
					getline(*in_file_, current_line_);
					break;
				}
			}
		}

		builder_->read_faces(zone_id_, first_index_, last_index_, bc_type_, face_type_, in_file_);
		set_operation_id(read_operation_id);
		break;

	case read_zone_39:
		zone_id_ = private_parse_number(10);
		private_parse_name();
		zone_name_ = private_parse_name();
		builder_->read_zones(zone_id_, zone_name_, in_file_);
		set_operation_id(read_operation_id);
		break;

	case read_zone_45:

		zone_id_ = private_parse_number(10);
		private_parse_name();
		zone_name_ = private_parse_name();
		builder_->read_zones(zone_id_, zone_name_, in_file_);
		set_operation_id(read_operation_id);
		break;

	case read_cell_tree:

		first_index_ = private_parse_number(16);
		last_index_ = private_parse_number(16);
		parent_id_ = private_parse_number(16);
		children_id_ = private_parse_number(16);
		for (int i = current_line_.length(); i > 0; --i)
		{
			char ch = current_line_[i];
			if (ch == ' ' || ch == '\0')
				continue;
			else if (ch == '(')
				break;
			else
			{
				getline(*in_file_, current_line_);
				break;
			}
		}
		builder_->read_cell_tree(first_index_, last_index_, parent_id_, children_id_, in_file_);
		set_operation_id(read_operation_id);
		break;


	case read_face_tree:
		first_index_ = private_parse_number(16);
		last_index_ = private_parse_number(16);
		parent_id_ = private_parse_number(16);
		children_id_ = private_parse_number(16);
		for (int i = current_line_.length(); i > 0; --i)
		{
			char ch = current_line_[i];
			if (ch == ' ' || ch == '\0')
				continue;
			else if (ch == '(')
				break;
			else
			{
				getline(*in_file_, current_line_);
				break;
			}
		}
		builder_->read_face_tree(first_index_, last_index_, parent_id_, children_id_, in_file_);
		set_operation_id(read_operation_id);
		break;
		break;
	default:

		set_operation_id(read_operation_id);
		break;
	}
	return;
}

psstInterfaceInputGrid * psstAnsysDirector::construct_grid(std::ifstream & in_stream)
{
	builder_->reset();
	in_file_ = &in_stream;
	while (!in_file_->eof())
	{
		build_step();
	}
	return builder_->get_grid();
}

size_t psstAnsysDirector::private_parse_number(size_t base)
{
	size_t destination = 0;
	size_t str_size = current_line_.size();

	char temp_ch;
	std::string buffer;
	bool is_empty_buffer = true;

	while (line_pos_ < str_size)
	{
		temp_ch = current_line_[line_pos_];
		if (base == 16)
		{
			if ((temp_ch >= '0' && temp_ch <= '9') || (temp_ch >= 'a' && temp_ch <= 'f') || (temp_ch >= 'A' && temp_ch <= 'F'))
			{
				buffer.push_back(temp_ch);
				is_empty_buffer = false;
			}
			else if (temp_ch == ' ' || temp_ch == '(' || temp_ch == ')')
			{
				if (!is_empty_buffer)
				{
					destination = std::stoi(buffer, 0, 16);
					break;
				}
			}
			else
			{
				break;
			}
		}
		else if (base == 10)
		{
			if (temp_ch >= '0' && temp_ch <= '9')
			{
				buffer.push_back(temp_ch);
				is_empty_buffer = false;
			}
			else if (temp_ch == ' ' || temp_ch == '(' || temp_ch == ')')
			{
				if (!is_empty_buffer)
				{
					destination = std::stoi(buffer, 0, 10);
					break;
				}
			}
			else
			{
				break;
			}
		}
		line_pos_++;
	}
	return destination;
}

std::string psstAnsysDirector::private_parse_name()
{
	size_t str_size = current_line_.size();

	char temp_ch;
	std::string buffer;
	bool is_empty_buffer = true;

	while (line_pos_ < str_size)
	{
		temp_ch = current_line_[line_pos_];
		if (temp_ch != ' ' && temp_ch != '(' && temp_ch != ')')
		{
			buffer.push_back(temp_ch);
			is_empty_buffer = false;
		}
		else if (!is_empty_buffer)
		{
			break;
		}
		line_pos_++;
	}
	return buffer;
}
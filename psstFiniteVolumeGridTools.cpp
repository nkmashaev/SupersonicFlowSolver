#include <exception>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include "psstVariable.h"
#include "psstFiniteVolumeGridTools.h"


/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceVertex																																													|
|	psstVertex																																															|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstVertex::psstVertex() :
	psstInterfaceVertex(),
	components_(),
	dimension_(0),
	id_(0)
{
	private_set_dimension(2);

	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
		components_[i] = 0.0;
}

psstVertex::psstVertex(size_t dimension) :
	psstInterfaceVertex(),
	components_(0),
	dimension_(0),
	id_(0)
{
	private_set_dimension(dimension);

	components_.resize(psstVertex::dimension());
	for (size_t i = 0; i < psstVertex::dimension(); ++i)
		components_[i] = 0.0;
}


psstVertex::psstVertex(const psstVertex & vertex) :
	psstInterfaceVertex(vertex),
	components_(0),
	dimension_(0),
	id_(0)
{
	private_set_dimension(vertex.dimension());

	components_.resize(psstVertex::dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = vertex.component(i);
	}
	id() = vertex.id();
}

psstVertex::psstVertex(const std::vector<double> & vertex) :
	psstInterfaceVertex(),
	components_(0),
	dimension_(0),
	id_(0)
{
	private_set_dimension(vertex.size());

	components_.resize(psstVertex::dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = vertex.at(i);
	}
}

psstVertex::psstVertex(const double * components, size_t dimension) :
	psstInterfaceVertex(),
	components_(0),
	dimension_(0),
	id_(0)
{
	if (components == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");

	private_set_dimension(dimension);
	components_.resize(psstVertex::dimension());
	for (size_t i = 0; i < psstVertex::dimension(); ++i)
	{
		components_[i] = components[i];
	}
}

void psstVertex::initialize()
{
	for (size_t i = 0; i < dimension(); ++i)
		components_[i] = 0.0;
}

psstInterfaceVertex * psstVertex::clone() const
{
	psstVertex * new_vertex = new psstVertex(*this);
	return new_vertex;
}

void psstVertex::copy(const psstInterfaceVertex & vertex)
{
	if (&vertex == this)
		return;

	private_set_dimension(vertex.dimension());
	
	components_.resize(psstVertex::dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = vertex.component(i);
	}
	id() = vertex.id();
}



const double & psstVertex::component(size_t i) const
{
	if (i >= dimension())
		throw std::out_of_range("Error: Out of range!");
	return components_[i];
}

psstVertex & psstVertex::operator=(const psstVertex & vertex)
{
	copy(vertex);
	return *this;
}

psstVertex::~psstVertex()
{}

void psstVertex::private_set_dimension(size_t dim)
{
	if (dim > 3)
		throw std::length_error("Error: Invalid dimension!");
	dimension_ = dim;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceVector																																													|
|	psstGeometryVector																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstGeometryVector::psstGeometryVector() :
	psstInterfaceVector(),
	components_(0),
	length_(0.0),
	dimension_(0)
{
	private_set_dimension(2);
	private_set_length(0.0);
	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
		components_[i] = 0.0;
}

psstGeometryVector::psstGeometryVector(size_t dim) :
	psstInterfaceVector(),
	components_(0),
	length_(0.0),
	dimension_(0)
{
	private_set_dimension(dim);
	private_set_length(0.0);
	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
		components_[i] = 0.0;
}

psstGeometryVector::psstGeometryVector(const psstGeometryVector & vector) :
	psstInterfaceVector(vector),
	components_(0),
	length_(0.0),
	dimension_(0)
{
	private_set_dimension(vector.dimension());
	private_set_length(vector.length());

	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = vector.component(i);
	}
}

psstGeometryVector::psstGeometryVector(const std::vector<double> & components) :
	psstInterfaceVector(),
	components_(0),
	length_(0.0),
	dimension_(0)
{
	private_set_dimension(components.size());
	
	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = components.at(i);
	}

	initialize();
}

psstGeometryVector::psstGeometryVector(const double * components, size_t dim) :
	psstInterfaceVector(),
	components_(0),
	length_(0.0),
	dimension_(0)
{
	if (components == 0)
		throw std::invalid_argument("Error: Expected data bun null_ptr found!");

	private_set_dimension(dim);

	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = components[i];
	}

	initialize();
}

psstGeometryVector::psstGeometryVector(const psstInterfaceVertex & start, const psstInterfaceVertex & end) :
	psstInterfaceVector(),
	components_(0),
	length_(0.0),
	dimension_(0)
{
	if (start.dimension() != end.dimension())
		throw std::domain_error("Error: Expected points with the same dimension!");
	private_set_dimension(start.dimension());

	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = end.component(i) - start.component(i);
	}

	initialize();
}

void psstGeometryVector::initialize()
{
	length_ = 0.0;
	for (size_t i = 0; i < dimension_; ++i)
	{
		length_ += component(i) * component(i);
	}
	length_ = std::sqrt(length_);
}

psstInterfaceVector * psstGeometryVector::clone() const
{
	psstGeometryVector * new_vector = new psstGeometryVector(*this);
	return new_vector;
}

void psstGeometryVector::copy(const psstInterfaceVector & vector)
{
	if (&vector == this)
		return;

	private_set_dimension(vector.dimension());
	private_set_length(vector.length());
	
	components_.resize(dimension());
	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] = vector.component(i);
	}
}

void psstGeometryVector::normalize()
{
	double eps = 1.0e-10;
	if (std::abs(length()) < eps)
		throw std::overflow_error("Error: Division by zero detected!");

	for (size_t i = 0; i < dimension(); ++i)
	{
		components_[i] /= length();
	}
	private_set_length(1.0);
}

double psstGeometryVector::dot_product(const psstInterfaceVector & vector) const
{
	if (vector.dimension() != dimension())
		throw std::domain_error("Error: Expected vectors with the same dimension!");

	double res = 0.0;
	for (size_t i = 0; i < dimension(); ++i)
		res += vector.component(i) * component(i);
	return res;
}

const double & psstGeometryVector::component(size_t i) const
{
	if (i >= dimension())
		throw std::out_of_range("Error: Out of range!");
	return components_[i];
}

psstGeometryVector::~psstGeometryVector()
{}

void psstGeometryVector::private_set_length(double length)
{
	if (length < 0.0)
		throw std::invalid_argument("Error: Expected positive length value!");
	length_ = length;
}

void psstGeometryVector::private_set_dimension(size_t dim)
{
	if (dim > 3)
		throw std::length_error("Error: Invalid dimension!");
	dimension_ = dim;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceFace																																													|
|	psstFace																																															|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstFace::psstFace() :
	psstInterfaceFace(),
	boundary_status_(false),
	square_(0.0),
	vertiñes_(),
	center_(),
	normal_vector_(),
	volumes_(),
	children_(),
	parent_(0),
	type_(0),
	id_(0)
{
	volumes_.resize(2);
	center_ = std::make_shared<psstVertex>();
	normal_vector_ = std::make_shared<psstGeometryVector>();
}

psstFace::psstFace(size_t dim) :
	psstInterfaceFace(),
	boundary_status_(false),
	square_(0.0),
	vertiñes_(),
	center_(),
	normal_vector_(),
	volumes_{ 0, 0 },
	children_(),
	parent_(0),
	id_(0)
{
	volumes_.resize(2);
	center_ = std::make_shared<psstVertex>(dim);
	normal_vector_ = std::make_shared<psstGeometryVector>(dim);
}

psstFace::psstFace(const psstFace & face) :
	psstInterfaceFace(face),
	boundary_status_(false),
	square_(0.0),
	vertiñes_(),
	center_(),
	normal_vector_(),
	volumes_{ 0, 0 },
	children_(),
	parent_(0),
	type_(0),
	id_(0)
{
	center_ = std::shared_ptr<psstInterfaceVertex>(face.center_->clone());
	normal_vector_ = std::shared_ptr<psstInterfaceVector>(face.normal_vector_->clone());

	boundary_status_ = face.is_boundary();
	square_ = face.square();
	vertiñes_ = face.vertiñes_;
	children_ = face.children_;
	parent_ = face.parent_;
	volumes_.resize(2);
	type_ = face.type();
	for (size_t i = 0; i < 2; ++i)
		volumes_[i] = face.neighbour_ptr(i);
	id() = face.id();
}

void psstFace::initialize()
{
	square_ = 0.0;
	for (size_t i = 0; i < dimension(); ++i)
	{
		normal_vector_->component(i) = 0.0;
		center_->component(i) = 0.0;
	}

	size_t dim = dimension();
	size_t vert_number = number_of_vertices();
	if (dim == 2)
	{
		if (vert_number != 2)
			throw std::logic_error("Error: Expected to vertices but extra data found!");
		double component = 0.0;
		
		for (size_t i = 0; i < dim; ++i)
		{
			center_->component(i) = 0.5 * (vertex(1)->component(i) + vertex(0)->component(i));
			component = vertex(1)->component(i) - vertex(0)->component(i);
			square_ += component * component;
		}
		square_ = std::sqrt(square_);
		
		
		normal_vector_->component(0) = vertex(0)->component(1) - vertex(1)->component(1);
		normal_vector_->component(1) = vertex(1)->component(0) - vertex(0)->component(0);
	}
	else if (dim == 3)
	{
		if (vert_number < 3)
			throw std::logic_error("Error: Expected at least 3 vertices!");

		double average[3] = {0.0, 0.0, 0.0};
		double average3[3] = { 0.0, 0.0, 0.0 };
		for (size_t i = 0; i < dim; ++i)
		{
			for (size_t j = 0; j < vert_number; ++j)
			{
				average[i] += vertex(j)->component(i);
			}
			average[i] /= vert_number;
			average3[i] = average[i] * 3;
		}

		//face center calculation
		double a[3] = { 0.0, 0.0, 0.0 };
		double b[3] = { 0.0, 0.0, 0.0 };
		double c[3] = { 0.0, 0.0, 0.0 };
		double temp_square = 0.0;
		size_t last_index = vert_number - 1;
		for (size_t i = 0; i < last_index; ++i)
		{
			
			for (size_t j = 0; j < dim; ++j)
			{
				a[j] = vertex(i)->component(j) - average[j];
				b[j] = vertex(i + 1)->component(j) - average[j];
			}
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];

			temp_square = 0.0;
			for (size_t j = 0; j < dim; ++j)
			{
				temp_square += c[j] * c[j];
			}
			temp_square = std::sqrt(temp_square);
			square_ += temp_square;

			for (size_t j = 0; j < dim; ++j)
			{
				center_->component(j) += temp_square * (a[j] + b[j] + average3[j]);
			}
		}
		
		for (size_t j = 0; j < dim; ++j)
		{
			a[j] = vertex(last_index)->component(j) - average[j];
			b[j] = vertex(0)->component(j) - average[j];
		}

		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];

		temp_square = 0.0;
		for (size_t j = 0; j < dim; ++j)
		{
			temp_square += c[j] * c[j];
		}
		temp_square = std::sqrt(temp_square);
		square_ += temp_square;

		double square3 = square_ * 3;
		for (size_t j = 0; j < dim; ++j)
		{
			center_->component(j) += temp_square * (a[j] + b[j] + average3[j]);
			center_->component(j) /= square3;
			normal_vector_->component(j) = c[j];
		}	
		square_ *= 0.5;
	}

	normal_vector_->initialize();
	normal_vector_->normalize();
}

psstInterfaceFace * psstFace::clone() const
{
	psstFace * new_face = new psstFace(*this);
	return new_face;
}

void psstFace::copy(const psstInterfaceFace & face)
{
	if (&face == this)
		return;
	center_ = std::shared_ptr<psstInterfaceVertex>(face.center_vertex()->clone());
	normal_vector_ = std::shared_ptr<psstInterfaceVector>(face.normal_vector()->clone());

	boundary_status_ = face.is_boundary();
	square_ = face.square();

	for (size_t i = 0; i < number_of_vertices(); ++i)
	{
		vertiñes_[i] = face.vertex_ptr(i);
	}

	for (size_t i = 0; i < 2; ++i)
		volumes_[i] = face.neighbour_ptr(i);

	parent_ = face.parent_ptr();
	for (size_t i = 0; i < face.number_of_children(); ++i)
	{
		children_[i] = face.child_ptr(i);
	}
	id() = face.id();
	set_type(face.type());
}

void psstFace::add_vertex(std::shared_ptr<psstInterfaceVertex> vertex)
{
	if (vertex.get() == 0)
		throw std::invalid_argument("Error: Expected vertex data but null_ptr found!");
	if (vertex->dimension() != dimension())
		throw std::invalid_argument("Error: Wrond dimension!");
	vertiñes_.push_back(vertex);
}

std::shared_ptr<psstInterfaceVertex> & psstFace::vertex(size_t i)
{
	if (i >= number_of_vertices())
		throw std::out_of_range("Error: Out of range!");
	return vertiñes_[i];
}

std::shared_ptr<psstInterfaceVertex> psstFace::vertex_ptr(size_t i) const
{
	if (i >= number_of_vertices())
		throw std::out_of_range("Error: Out of range!");
	return vertiñes_[i];
}

std::shared_ptr<const psstInterfaceVertex> psstFace::vertex(size_t i) const
{
	if (i >= number_of_vertices())
		throw std::out_of_range("Error: Out of range!");
	return vertiñes_[i];
}

size_t psstFace::number_of_vertices() const
{
	return vertiñes_.size();
}

std::shared_ptr<psstInterfaceVolume> & psstFace::neighbour(size_t i)
{
	if (i >= 2)
		throw std::out_of_range("Error: Out of range!");
	return volumes_[i];
}

std::shared_ptr<psstInterfaceVolume> psstFace::neighbour_ptr(size_t i) const
{
	if (i >= 2)
		throw std::out_of_range("Error: Out of range!");
	return volumes_[i];
}

std::shared_ptr<const psstInterfaceVolume> psstFace::neighbour(size_t i) const
{
	if (i >= 2)
		throw std::out_of_range("Error: Out of range!");
	return volumes_[i];
}

std::shared_ptr<psstInterfaceVolume> & psstFace::neighbour_for(std::shared_ptr<psstInterfaceVolume> vol)
{
	if (vol.get() != volumes_[0].get() && vol.get() != volumes_[1].get())
		throw std::invalid_argument("Error: Input volume expected to be face neighbour!");
	if (vol.get() == volumes_[0].get())
		return volumes_[1];
	return volumes_[0];
}

std::shared_ptr<const psstInterfaceVolume> psstFace::neighbour_for(std::shared_ptr<const psstInterfaceVolume> vol) const
{
	if (vol.get() != volumes_[0].get() && vol.get() != volumes_[1].get())
		throw std::invalid_argument("Error: Input volume expected to be face neighbour!");
	if (vol.get() == volumes_[0].get())
		return volumes_[1];
	return volumes_[0];
}

std::shared_ptr<psstInterfaceVolume> & psstFace::neighbour_for(psstInterfaceVolume * vol)
{
	if (vol != volumes_[0].get() && vol != volumes_[1].get())
		throw std::invalid_argument("Error: Input volume expected to be face neighbour!");
	if (vol == volumes_[0].get())
		return volumes_[1];
	return volumes_[0];
}

std::shared_ptr<const psstInterfaceVolume> psstFace::neighbour_for(const psstInterfaceVolume * vol) const
{
	if (vol != volumes_[0].get() && vol != volumes_[1].get())
		throw std::invalid_argument("Error: Input volume expected to be face neighbour!");
	if (vol == volumes_[0].get())
		return volumes_[1];
	return volumes_[0];
}

void psstFace::add_child(std::shared_ptr<psstInterfaceFace> child)
{
	if (child.get() == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");
	children_.push_back(child);
}

psstFace & psstFace::operator=(const psstFace & face)
{
	copy(face);
	return *this;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceVolume																																													|
|	psstVolume																																															|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstVolume::psstVolume() :
	boundary_status_(false),
	volume_(0.0),
	neighbours_(),
	type_(0),
	refine_level_(0),
	center_(0),
	variable_(0),
	children_(),
	parent_(0),
	id_(0)
{
	center_ = std::make_shared<psstVertex>();
	variable_ = std::make_shared<psstVariable>();
}

psstVolume::psstVolume(size_t dim) :
	boundary_status_(false),
	volume_(0.0),
	neighbours_(),
	type_(0),
	refine_level_(0),
	center_(0),
	variable_(0),
	children_(),
	parent_(0),
	id_(0)
{
	center_ = std::make_shared<psstVertex>(dim);
	variable_ = std::make_shared<psstVariable>(dim);
}

psstVolume::psstVolume(const psstVolume & vol) :
	boundary_status_(false),
	volume_(0.0),
	neighbours_(),
	type_(0),
	refine_level_(0),
	center_(0),
	variable_(0),
	children_(),
	parent_(0),
	id_(0)
{
	center_ = std::shared_ptr<psstInterfaceVertex>(vol.center_->clone());
	variable_ = std::shared_ptr<psstInterfaceVariable>(vol.variable_->clone());
	neighbours_ = vol.neighbours_;

	boundary_status_ = vol.is_boundary();
	type_ = vol.type();
	refine_level_ = vol.refine_level();

	children_ = vol.children_;
	parent_ = vol.parent_;
	id_ = vol.id_;
}

void psstVolume::initialize()
{
	volume_ = 0.0;
	size_t dim = dimension();
	size_t face_numb = number_of_faces();

	for (size_t i = 0; i < dim; ++i)
	{
		center_->component(i) = 0.0;
	}
	
	for (size_t i = 0; i < face_numb; ++i)
	{
		neighbours_[i].volume_neighbour = neighbours_[i].face_neighbour->neighbour_for(this);
	}

	if (is_boundary())
	{
		center_ = std::shared_ptr<psstInterfaceVertex>(neighbours_[0].face_neighbour->center_vertex()->clone());
		return;
	}


	if (dim == 2)
	{
		double average[2] = { 0.0, 0.0 };
		double average3[2] = { 0.0, 0.0 };
		for (size_t j = 0; j < dim; ++j)
		{
			for (size_t i = 0; i < face_numb; ++i)
			{
				average[j] += face(i)->center_vertex()->component(j);
			}
			average[j] /= face_numb;
			average3[j] = average[j] * 3;
		}

		double temp_volume = 0.0;
		double a[2] = { 0.0, 0.0 };
		double b[2] = { 0.0, 0.0 };
		for (size_t i = 0; i < face_numb; ++i)
		{
			for (size_t k = 0; k < dim; ++k)
			{
				a[k] = face(i)->vertex(0)->component(k) - average[k];
				b[k] = face(i)->vertex(1)->component(k) - average[k];
			}
			temp_volume = std::abs(a[0] * b[1] - a[1] * b[0]);
			volume_ += temp_volume;
			for (size_t k = 0; k < dim; ++k)
			{
				center_->component(k) += temp_volume * (a[k] + b[k] + average3[k]);
			}
		}

		double volume3 = volume_ * 3;
		for (size_t k = 0; k < dim; ++k)
		{
			center_->component(k) /= volume3;
		}
		volume_ /= 2;

	}
	else if (dim == 3)
	{
		double average[3] = { 0.0, 0.0, 0.0 };
		double average4[3] = { 0.0, 0.0, 0.0 };
		for (size_t j = 0; j < dim; ++j)
		{
			for (size_t i = 0; i < face_numb; ++i)
			{
				average[j] += face(i)->center_vertex()->component(j);
			}
			average[j] /= face_numb;
			average4[j] = average[j] * 4;
		}

		double temp_volume = 0.0;
		double a[3] = { 0.0, 0.0, 0.0 };
		double b[3] = { 0.0, 0.0, 0.0 };
		double c[3] = { 0.0, 0.0, 0.0 };
		size_t next = 0;
		size_t vert_numb = 0;
		std::shared_ptr<const psstInterfaceVertex> center = 0;
		for (size_t i = 0; i < face_numb; ++i)
		{
			vert_numb = face(i)->number_of_vertices();
			for (size_t j = 0; j < vert_numb; ++j)
			{
				next = j + 1;
				if (next == vert_numb)
					next = 0;

				center = face(i)->center_vertex();
				for (size_t k = 0; k < dim; ++k)
				{
					a[k] = face(i)->vertex(j)->component(k) - average[k];
					b[k] = face(i)->vertex(next)->component(k) - average[k];
					c[k] = center->component(k) - average[k];
				}

				temp_volume = std::abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]));
				volume_ += temp_volume;
				for (size_t k = 0; k < dimension(); ++k)
				{
					center_->component(k) += temp_volume * (a[k] + b[k] + c[k] + average4[k]);
				}
			}
		}	

		double volume4 = volume_ * 4.0;
		for (size_t k = 0; k < dimension(); ++k)
		{
			center_->component(k) /= volume4;
		}
		volume_ /= 6.0;
	}
	
}

psstInterfaceVolume * psstVolume::clone() const
{
	psstVolume * new_vol = new psstVolume(*this);
	return new_vol;
}

void psstVolume::copy(const psstInterfaceVolume & vol)
{
	if (&vol == this)
		return;

	center_ = std::shared_ptr<psstInterfaceVertex>(vol.center_vertex()->clone());
	variable_ = std::shared_ptr<psstInterfaceVariable>(vol.variable()->clone());

	for (size_t i = 0; i < number_of_faces(); ++i)
	{
		face(i) = vol.face_ptr(i);
		neighbour(i) = vol.neighbour_ptr(i);
	}

	boundary_status_ = vol.is_boundary();
	type_ = vol.type();
	refine_level_ = vol.refine_level();

	parent_ = vol.parent_ptr();
	for (size_t i = 0; i < vol.number_of_children(); ++i)
	{
		children_[i] = vol.child_ptr(i);
	}
	id() = vol.id();
}

void psstVolume::add_face(std::shared_ptr<psstInterfaceFace> face)
{
	if (face.get() == 0)
		throw std::invalid_argument("Error: Expected face data but null_ptr found!");
	psstNeighbour neighbour;
	neighbour.face_neighbour = face;
	neighbour.volume_neighbour = 0;
	neighbours_.push_back(neighbour);
}

std::shared_ptr<psstInterfaceFace> & psstVolume::face(size_t i)
{
	return neighbours_.at(i).face_neighbour;
}

std::shared_ptr<psstInterfaceFace> psstVolume::face_ptr(size_t i) const
{
	return neighbours_.at(i).face_neighbour;
}

std::shared_ptr<const psstInterfaceFace> psstVolume::face(size_t i) const
{
	return neighbours_.at(i).face_neighbour;
}

std::shared_ptr<psstInterfaceVolume> & psstVolume::neighbour(size_t i)
{
	return neighbours_.at(i).volume_neighbour;
}

std::shared_ptr<psstInterfaceVolume> psstVolume::neighbour_ptr(size_t i) const
{
	return neighbours_.at(i).volume_neighbour;
}

std::shared_ptr<const psstInterfaceVolume> psstVolume::neighbour(size_t i) const
{
	return neighbours_.at(i).volume_neighbour;
}

void psstVolume::add_child(std::shared_ptr<psstInterfaceVolume> child)
{
	if (child.get() == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");
	children_.push_back(child);
}

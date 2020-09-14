#ifndef _PSST_FINITE_VOLUME_GRID_TOOLS_H_
#define _PSST_FINITE_VOLUME_GRID_TOOLS_H_

#include<vector>
#include<memory>
class psstInterfaceVertex
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceVertex * clone() const = 0;
	virtual void copy(const psstInterfaceVertex &) = 0;

	virtual size_t & id() = 0;
	virtual const size_t & id() const = 0;

	virtual size_t dimension() const = 0;
	virtual double & component(size_t i) = 0;
	virtual const double & component(size_t i) const = 0;

	virtual ~psstInterfaceVertex() {}
};

class psstVertex : public psstInterfaceVertex
{
private:
	std::vector<double> components_;
	size_t dimension_;
	size_t id_;
	void private_set_dimension(size_t);
public:
	psstVertex();
	psstVertex(size_t dimension);
	psstVertex(const psstVertex &);
	psstVertex(const std::vector<double> &);
	psstVertex(const double * components, size_t dimenstion);

	virtual void initialize();
	virtual psstInterfaceVertex * clone() const;
	virtual void copy(const psstInterfaceVertex &);

	virtual size_t & id() { return id_; };
	virtual const size_t & id() const { return id_; };

	virtual size_t dimension() const { return dimension_; };
	virtual double & component(size_t i) { return const_cast<double &>( static_cast<const psstVertex &>(*this).component(i)); };
	virtual const double & component(size_t i) const;

	psstVertex & operator=(const psstVertex &);
	virtual ~psstVertex();
};

class psstInterfaceVector
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceVector* clone() const = 0;
	virtual void copy(const psstInterfaceVector &) = 0;

	virtual size_t dimension() const = 0;
	virtual void normalize() = 0;
	virtual double length() const = 0;
	virtual double dot_product(const psstInterfaceVector &) const = 0;
	virtual double & component(size_t i) = 0;
	virtual const double & component(size_t i) const = 0;
	virtual ~psstInterfaceVector() {};
};

class psstGeometryVector : public psstInterfaceVector
{
private:
	std::vector<double> components_;
	double length_;
	size_t dimension_;

	void private_set_dimension(size_t);
	virtual void private_set_length(double length);
public:
	psstGeometryVector();
	psstGeometryVector(size_t dimension);
	psstGeometryVector(const psstGeometryVector &);
	psstGeometryVector(const std::vector<double> &);
	psstGeometryVector(const double *, size_t);
	psstGeometryVector(const psstInterfaceVertex & start, const psstInterfaceVertex & end);

	virtual void initialize();
	virtual psstInterfaceVector * clone() const;
	virtual void copy(const psstInterfaceVector &);

	virtual size_t dimension() const { return dimension_; };
	virtual void normalize();
	virtual double length() const { return length_; };
	virtual double dot_product(const psstInterfaceVector &) const;
	virtual double & component(size_t i) { return const_cast<double &>(static_cast<const psstGeometryVector &>(*this).component(i)); };
	virtual const double & component(size_t i) const;

	virtual ~psstGeometryVector();
};

class psstInterfaceVolume;
class psstInterfaceFace
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceFace * clone() const = 0;
	virtual void copy(const psstInterfaceFace &) = 0;

	virtual size_t dimension() const = 0;
	virtual size_t & id() = 0;
	virtual const size_t & id() const = 0;

	virtual void add_vertex(std::shared_ptr<psstInterfaceVertex>) = 0;
	virtual std::shared_ptr<psstInterfaceVertex> & vertex(size_t i) = 0;
	virtual std::shared_ptr<psstInterfaceVertex> vertex_ptr(size_t i) const = 0;
	virtual std::shared_ptr<const psstInterfaceVertex> vertex(size_t i) const = 0;
	virtual size_t number_of_vertices() const = 0;

	virtual std::shared_ptr<psstInterfaceVolume> & neighbour(size_t i)= 0;
	virtual std::shared_ptr<psstInterfaceVolume> neighbour_ptr(size_t i) const = 0;
	virtual std::shared_ptr<const psstInterfaceVolume> neighbour(size_t) const = 0;
	virtual std::shared_ptr<psstInterfaceVolume> & neighbour_for(std::shared_ptr<psstInterfaceVolume>) = 0;
	virtual std::shared_ptr<const psstInterfaceVolume> neighbour_for(std::shared_ptr<const psstInterfaceVolume>) const = 0;
	virtual std::shared_ptr<psstInterfaceVolume> & neighbour_for(psstInterfaceVolume *) = 0;
	virtual std::shared_ptr<const psstInterfaceVolume> neighbour_for(const psstInterfaceVolume *) const = 0;
	virtual double square() const = 0;

	virtual std::shared_ptr<const psstInterfaceVertex> center_vertex() const = 0;
	virtual std::shared_ptr<const psstInterfaceVector> normal_vector() const = 0;
	virtual bool is_boundary() const = 0;
	virtual void set_boundary_status(bool status) = 0;

	virtual size_t type() const = 0;
	virtual void set_type(size_t type) = 0;

	virtual size_t number_of_children() const = 0;
	virtual void add_child(std::shared_ptr<psstInterfaceFace>) = 0;
	virtual std::shared_ptr<psstInterfaceFace> & child(size_t i) = 0;
	virtual std::shared_ptr<psstInterfaceFace> child_ptr(size_t i) const = 0;
	virtual std::shared_ptr<const psstInterfaceFace> child(size_t i) const = 0;

	virtual bool have_a_parent() const = 0;
	virtual std::shared_ptr<psstInterfaceFace> & parent() = 0;
	virtual std::shared_ptr<psstInterfaceFace> parent_ptr() const = 0 ;
	virtual std::shared_ptr<const psstInterfaceFace> parent() const = 0;

	virtual ~psstInterfaceFace() {};
};

class psstFace : public psstInterfaceFace
{
private:
	bool boundary_status_;
	double square_;
	size_t type_;
	size_t id_;
	std::vector<std::shared_ptr<psstInterfaceVertex>> vertiñes_;
	std::shared_ptr<psstInterfaceVertex> center_;
	std::shared_ptr<psstInterfaceVector> normal_vector_;
	std::vector<std::shared_ptr<psstInterfaceVolume>> volumes_; //0 - right volume, 1 - left volume

	std::vector<std::shared_ptr<psstInterfaceFace>> children_;
	std::shared_ptr<psstInterfaceFace> parent_;
public:
	psstFace();
	psstFace(size_t dim);
	psstFace(const psstFace &);

	virtual void initialize();
	virtual psstInterfaceFace * clone() const;
	virtual void copy(const psstInterfaceFace &);

	virtual size_t dimension() const { return center_->dimension(); };
	virtual size_t & id() { return id_; };
	virtual const size_t & id() const { return id_; };

	virtual void add_vertex(std::shared_ptr<psstInterfaceVertex>);
	virtual std::shared_ptr<psstInterfaceVertex> & vertex(size_t i);
	virtual std::shared_ptr<psstInterfaceVertex> vertex_ptr(size_t i) const;
	virtual std::shared_ptr<const psstInterfaceVertex> vertex(size_t i) const;
	virtual size_t number_of_vertices() const;

	virtual std::shared_ptr<psstInterfaceVolume> & neighbour(size_t i);
	virtual std::shared_ptr<psstInterfaceVolume>  neighbour_ptr(size_t i) const;
	virtual std::shared_ptr<const psstInterfaceVolume>  neighbour(size_t) const;
	virtual std::shared_ptr<psstInterfaceVolume> & neighbour_for(std::shared_ptr<psstInterfaceVolume>);
	virtual std::shared_ptr<const psstInterfaceVolume>  neighbour_for(std::shared_ptr<const psstInterfaceVolume>) const;
	virtual std::shared_ptr<psstInterfaceVolume> & neighbour_for(psstInterfaceVolume *);
	virtual std::shared_ptr<const psstInterfaceVolume> neighbour_for(const psstInterfaceVolume *) const;
	virtual double square() const { return square_; };

	virtual std::shared_ptr<const psstInterfaceVertex> center_vertex() const { return center_; };
	virtual std::shared_ptr<const psstInterfaceVector> normal_vector() const { return normal_vector_; };
	virtual bool is_boundary() const { return boundary_status_; };
	virtual void set_boundary_status(bool status) { boundary_status_ = status; };

	virtual size_t type() const { return type_; };
	virtual void set_type(size_t type) { type_ = type; };

	virtual size_t number_of_children() const { return children_.size(); };
	virtual void add_child(std::shared_ptr<psstInterfaceFace>);
	virtual std::shared_ptr<psstInterfaceFace> & child(size_t i) { return children_.at(i); };
	virtual std::shared_ptr<psstInterfaceFace> child_ptr(size_t i) const { return children_.at(i); };
	virtual std::shared_ptr<const psstInterfaceFace> child(size_t i) const { return children_.at(i); };

	virtual bool have_a_parent() const { return parent_.get() != 0; };
	virtual std::shared_ptr<psstInterfaceFace> & parent() { return parent_; };
	virtual std::shared_ptr<psstInterfaceFace> parent_ptr() const { return parent_; };
	virtual std::shared_ptr<const psstInterfaceFace> parent() const { return parent_; };

	psstFace & operator=(const psstFace &);
	virtual ~psstFace() {};

};

class psstInterfaceVariable;
class psstInterfaceVolume
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceVolume * clone() const = 0;
	virtual void copy(const psstInterfaceVolume &) = 0;

	virtual size_t dimension() const = 0;
	virtual size_t & id() = 0;
	virtual const size_t & id() const = 0;

	virtual double volume() const = 0;

	virtual void add_face(std::shared_ptr<psstInterfaceFace>) = 0;
	virtual std::shared_ptr<psstInterfaceFace> & face(size_t i) = 0;
	virtual std::shared_ptr<psstInterfaceFace> face_ptr(size_t i) const = 0;
	virtual std::shared_ptr<const psstInterfaceFace> face(size_t i) const = 0;
	virtual std::shared_ptr<psstInterfaceVolume> & neighbour(size_t i) = 0;
	virtual std::shared_ptr<psstInterfaceVolume> neighbour_ptr(size_t i) const = 0;
	virtual std::shared_ptr<const psstInterfaceVolume> neighbour(size_t i) const = 0;
	virtual size_t number_of_faces() const = 0;

	virtual std::shared_ptr<psstInterfaceVariable> & variable() = 0;
	virtual std::shared_ptr<const psstInterfaceVariable> variable() const = 0;

	virtual std::shared_ptr<psstInterfaceVertex> center_vertex() const =0;
	virtual bool is_boundary() const = 0;
	virtual void set_boundary_status(bool status) = 0;

	virtual size_t type() const = 0;
	virtual void set_type(size_t type) = 0;
	virtual size_t refine_level() const = 0;
	virtual void set_refine_level(size_t lvl) = 0;

	virtual size_t number_of_children() const = 0;
	virtual void add_child(std::shared_ptr<psstInterfaceVolume>) = 0;
	virtual std::shared_ptr<psstInterfaceVolume> & child(size_t i) = 0;
	virtual std::shared_ptr<psstInterfaceVolume> child_ptr(size_t i) const = 0;
	virtual std::shared_ptr<const psstInterfaceVolume> child(size_t i) const = 0;

	virtual bool have_a_parent() const = 0;
	virtual std::shared_ptr<psstInterfaceVolume> & parent() = 0;
	virtual std::shared_ptr<psstInterfaceVolume> parent_ptr() const = 0;
	virtual std::shared_ptr<const psstInterfaceVolume> parent() const = 0;

	virtual ~psstInterfaceVolume() {};
};

struct psstNeighbour
{
public:
	std::shared_ptr<psstInterfaceVolume> volume_neighbour;
	std::shared_ptr<psstInterfaceFace>  face_neighbour;
};

class psstVolume : public psstInterfaceVolume
{
private:
	bool boundary_status_;
	double volume_;
	std::vector<psstNeighbour> neighbours_;
	size_t id_;
	size_t type_;
	size_t refine_level_;
	std::shared_ptr<psstInterfaceVertex> center_;
	std::shared_ptr<psstInterfaceVariable> variable_;

	std::vector<std::shared_ptr<psstInterfaceVolume>> children_;
	std::shared_ptr<psstInterfaceVolume> parent_;
public:
	psstVolume();
	psstVolume(size_t dim);
	psstVolume(const psstVolume &);

	virtual void initialize();
	virtual psstInterfaceVolume * clone() const;
	virtual void copy(const psstInterfaceVolume &);

	virtual size_t dimension() const { return center_->dimension(); };
	virtual size_t & id() { return id_; };
	virtual const size_t & id() const { return id_; };

	virtual double volume() const { return volume_; };

	virtual void add_face(std::shared_ptr<psstInterfaceFace>);
	
	virtual std::shared_ptr<psstInterfaceFace> & face(size_t i);
	virtual std::shared_ptr<psstInterfaceFace> face_ptr(size_t i) const;
	virtual std::shared_ptr<const psstInterfaceFace> face(size_t i) const;
	virtual std::shared_ptr<psstInterfaceVolume> & neighbour(size_t i);
	virtual std::shared_ptr<psstInterfaceVolume> neighbour_ptr(size_t i) const;
	virtual std::shared_ptr<const psstInterfaceVolume> neighbour(size_t i) const;
	virtual size_t number_of_faces() const { return neighbours_.size(); };

	virtual std::shared_ptr<psstInterfaceVariable> & variable() { return variable_; };
	virtual std::shared_ptr<const psstInterfaceVariable> variable() const { return variable_; };

	virtual std::shared_ptr<psstInterfaceVertex> center_vertex() const { return center_; };
	virtual bool is_boundary() const { return boundary_status_; };
	virtual void set_boundary_status(bool status) { boundary_status_ = status; };

	virtual size_t type() const { return type_; };
	virtual void set_type(size_t type) { type_ = type; };
	virtual size_t refine_level() const { return refine_level_; };
	virtual void set_refine_level(size_t lvl) { refine_level_ = lvl; };

	virtual size_t number_of_children() const { return children_.size(); };
	virtual void add_child(std::shared_ptr<psstInterfaceVolume>);
	virtual std::shared_ptr<psstInterfaceVolume> & child(size_t i) { return children_.at(i); };
	virtual std::shared_ptr<psstInterfaceVolume> child_ptr(size_t i) const { return children_.at(i); };
	virtual std::shared_ptr<const psstInterfaceVolume> child(size_t i) const { return children_.at(i); };

	virtual bool have_a_parent() const { return parent_.get() != 0; };
	virtual std::shared_ptr<psstInterfaceVolume> & parent() { return parent_; };
	virtual std::shared_ptr<psstInterfaceVolume> parent_ptr() const { return parent_; };
	virtual std::shared_ptr<const psstInterfaceVolume> parent() const { return parent_; };

	psstInterfaceVolume & operator=(const psstInterfaceVolume & vol) { copy(vol); return *this; };
	virtual ~psstVolume() {};
};

struct psstFaceStorage
{
public:
	std::vector<std::shared_ptr<psstInterfaceFace>> faces_;
	size_t first_index_;
	size_t last_index_;
	size_t bc_type_id_;
	size_t children_zone_id_;
	size_t parent_zone_id_;
	std::string name_;
};

struct psstBoundary;
struct psstVolumeStorage
{
public:
	std::vector<std::shared_ptr<psstInterfaceVolume>> volumes_;
	size_t first_index_;
	size_t last_index_;
	size_t bc_type_id_;
	size_t children_zone_id_;
	size_t parent_zone_id_;
	std::shared_ptr<psstBoundary> condition_;
	std::string name_;
};


#endif
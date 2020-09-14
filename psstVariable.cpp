#include "psstVariable.h"
#include "psstFiniteVolumeGridTools.h"
#include <exception>
#include <stdexcept>

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceVariable																																												|
|	psstVariable																																														|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstVariable::psstVariable() :
	psstInterfaceVariable(),
	velocity_(0),
	temperature_(0.0),
	grad_velocity_(0),
	grad_temperature_(0),
	grad_pressure_(0),
	grad_density_(0),
	pressure_(0.0),
	density_(0.0),
	time_step_(0.0),
	dimension_(0)
{
	private_set_dimension(2);
	velocity_.resize(dimension());
	size_t size = dimension() + 2;
	residual_.resize(size);
	w_.resize(size);
	prev_w_.resize(size);
	grad_velocity_.resize(dimension());

	for (size_t j = 0; j < dimension(); ++j)
	{
		grad_velocity_[j] = std::make_shared<psstGeometryVector>(dimension());
	}
	grad_temperature_ = std::make_shared<psstGeometryVector>();
	grad_pressure_ = std::make_shared<psstGeometryVector>();
	grad_density_ = std::make_shared<psstGeometryVector>();

	temperature_ = 0.0;
	pressure_ = 0.0;
	density_ = 0.0;
	time_step_ = 0.0;
	for (size_t j = 0; j < dimension(); ++j)
	{
		velocity_[j] = 0.0;
		grad_temperature_->component(j) = 0.0;
		grad_pressure_->component(j) = 0.0;
		grad_density_->component(j) = 0.0;
		for (size_t i = 0; i < dimension(); ++i)
			grad_velocity_[i]->component(j) = 0.0;
		residual_[j] = 0.0;
		w_[j] = 0.0;
	}

	for (size_t j = dimension(); j < size; ++j)
	{
		residual_[j] = 0.0;
		w_[j] = 0.0;
		prev_w_[j] = 0.0;
	}
}

psstVariable::psstVariable(size_t dim) :
	psstInterfaceVariable(),
	velocity_(0),
	temperature_(0.0),
	grad_velocity_(0),
	grad_temperature_(0),
	grad_pressure_(0),
	grad_density_(0),
	pressure_(0.0),
	density_(0.0),
	time_step_(0.0),
	dimension_(0)
{
	private_set_dimension(dim);
	velocity_.resize(dimension());
	size_t size = dimension() + 2;
	residual_.resize(size);
	w_.resize(size);
	prev_w_.resize(size);
	grad_velocity_.resize(dimension());

	for (size_t j = 0; j < dimension(); ++j)
	{
		grad_velocity_[j] = std::make_shared<psstGeometryVector>(dimension());
	}
	grad_temperature_ = std::make_shared<psstGeometryVector>(dimension());
	grad_pressure_ = std::make_shared<psstGeometryVector>(dimension());
	grad_density_ = std::make_shared<psstGeometryVector>(dimension());

	temperature_ = 0.0;
	pressure_ = 0.0;
	density_ = 0.0;
	time_step_ = 0.0;
	for (size_t j = 0; j < dimension(); ++j)
	{
		velocity_[j] = 0.0;
		grad_temperature_->component(j) = 0.0;
		grad_pressure_->component(j) = 0.0;
		grad_density_->component(j) = 0.0;
		for (size_t i = 0; i < dimension(); ++i)
			grad_velocity_[i]->component(j) = 0.0;
		residual_[j] = 0.0;
		w_[j] = 0.0;
	}

	for (size_t j = dimension(); j < size; ++j)
	{
		residual_[j] = 0.0;
		w_[j] = 0.0;
		prev_w_[j] = 0.0;
	}
}

psstVariable::psstVariable(const psstVariable & var) :
	psstInterfaceVariable(var),
	velocity_(),
	temperature_(0.0),
	grad_velocity_(0),
	grad_temperature_(0),
	grad_pressure_(0),
	grad_density_(0),
	pressure_(0.0),
	density_(0.0),
	time_step_(0.0),
	dimension_(0),
	w_(),
	residual_()
{
	size_t dim = var.dimension();
	private_set_dimension(dim);
	velocity_.resize(dim);
	size_t size = dim + 2;
	residual_.resize(size);
	w_.resize(size);
	prev_w_.resize(size);
	grad_velocity_.resize(dim);

	for (size_t j = 0; j < dim; ++j)
	{
		grad_velocity_[j] = std::make_shared<psstGeometryVector>(dimension());
	}
	grad_temperature_ = std::make_shared<psstGeometryVector>(dimension());
	grad_pressure_ = std::make_shared<psstGeometryVector>(dimension());
	grad_density_ = std::make_shared<psstGeometryVector>(dimension());

	temperature_ = var.temperature();
	pressure_ = var.pressure();
	density_ = var.density();
	time_step_ = var.time_step();
	for (size_t j = 0; j < dim; ++j)
	{
		velocity_[j] = var.velocity(j);
		grad_temperature_->component(j) = var.grad_temperature()->component(j);
		grad_pressure_->component(j) = var.grad_pressure()->component(j);
		grad_density_->component(j) = var.grad_density()->component(j);
		for (size_t i = 0; i < dim; ++i)
			grad_velocity_[i]->component(j) = var.grad_velocity(i)->component(j);
	}

	for (size_t j = 0; j < size; ++j)
	{
		residual_[j] = var.residual(j);
		w_[j] = var.w(j);
		prev_w_[j] = var.prev_w_[j];
	}
}

void psstVariable::initialize()
{
	temperature_ = 0.0;
	pressure_ = 0.0;
	density_ = 0.0;
	time_step_ = 0.0;

	for (size_t j = 0; j < dimension(); ++j)
	{
		velocity_[j] = 0.0;
		grad_temperature_->component(j) = 0.0;
		grad_pressure_->component(j) = 0.0;
		grad_density_->component(j) = 0.0;
		for (size_t i = 0; i < dimension(); ++i)
			grad_velocity_[i]->component(j) = 0.0;
		residual_[j] = 0.0;
		w_[j] = 0.0;
	}

	size_t size = dimension() + 2;
	for (size_t j = dimension(); j < size; ++j)
	{
		residual_[j] = 0.0;
		w_[j] = 0.0;
	}
}

psstInterfaceVariable * psstVariable::clone() const
{
	psstVariable * var = new psstVariable(*this);
	return var;
}

void psstVariable::copy(const psstInterfaceVariable & var)
{
	if (&var == this)
		return;
	private_set_dimension(var.dimension());

	size_t dim = dimension();
	size_t size = dimension() + 2;
	
	velocity_.resize(dim);
	residual_.resize(size);
	w_.resize(size);
	grad_velocity_.resize(dim);

	for (size_t j = 0; j < dim; ++j)
	{
		grad_velocity_[j] = std::make_shared<psstGeometryVector>(dim);
	}
	grad_temperature_ = std::make_shared<psstGeometryVector>();
	grad_pressure_ = std::make_shared<psstGeometryVector>();
	grad_density_ = std::make_shared<psstGeometryVector>();

	temperature_ = var.temperature();
	pressure_ = var.pressure();
	density_ = var.density();
	time_step_ = var.time_step();
	for (size_t j = 0; j < dim; ++j)
	{
		velocity_[j] = var.velocity(j);
		grad_temperature_->component(j) = var.grad_temperature()->component(j);
		grad_pressure_->component(j) = var.grad_pressure()->component(j);
		grad_density_->component(j) = var.grad_density()->component(j);
		for (size_t i = 0; i < dim; ++i)
			grad_velocity_[i]->component(j) = var.grad_velocity(i)->component(j);
		residual_[j] = var.residual(j);
		w_[j] = var.w(j);
	}

	for (size_t j = dim; j < size; ++j)
	{
		residual_[j] = var.residual(j);
		w_[j] = var.w(j);
		prev_w_[j] = var.prev_w(j);
	}
}

double & psstVariable::velocity(size_t i)
{
	return velocity_[i];
}

const double & psstVariable::velocity(size_t i) const
{
	return velocity_[i];
}

double & psstVariable::temperature()
{
	return temperature_;
}

const double & psstVariable::temperature() const
{
	return temperature_;
}

double & psstVariable::pressure()
{
	return pressure_;
}

const double & psstVariable::pressure() const
{
	return pressure_;
}

double & psstVariable::density()
{
	return density_;
}

const double & psstVariable::density() const
{
	return density_;
}

std::shared_ptr<psstInterfaceVector> & psstVariable::grad_velocity(size_t i)
{
	return grad_velocity_[i];
}

std::shared_ptr<const psstInterfaceVector> psstVariable::grad_velocity(size_t i) const
{
	return grad_velocity_[i];
}

std::shared_ptr<psstInterfaceVector> & psstVariable::grad_temperature()
{
	return grad_temperature_;
}

std::shared_ptr<const psstInterfaceVector> psstVariable::grad_temperature() const
{
	return grad_temperature_;
}


std::shared_ptr<psstInterfaceVector>&  psstVariable::grad_pressure()
{
	return grad_pressure_;
}

std::shared_ptr<const psstInterfaceVector> psstVariable::grad_pressure() const
{
	return grad_pressure_;
}

std::shared_ptr<psstInterfaceVector> & psstVariable::grad_density()
{
	return grad_density_;
}

std::shared_ptr<const psstInterfaceVector> psstVariable::grad_density() const
{
	return grad_density_;
}

void psstVariable::nullify_grad()
{
	size_t dim = dimension();
	for (size_t j = 0; j < dim; ++j)
	{
		grad_temperature_->component(j) = 0.0;
		grad_pressure_->component(j) = 0.0;
		grad_density_->component(j) = 0.0;
		for (size_t i = 0; i < dim; ++i)
			grad_velocity_[i]->component(j) = 0.0;
	}
	grad_temperature_->initialize();
	grad_pressure_->initialize();
	grad_density_->initialize();
	for (size_t i = 0; i < dim; ++i)
		grad_velocity_[i]->initialize();
}

double & psstVariable::w(size_t i)
{
	return w_[i];
}

const double & psstVariable::w(size_t i) const
{
	return w_[i];
}

double & psstVariable::residual(size_t i)
{
	return residual_[i];
}

const double & psstVariable::residual(size_t i) const
{
	return residual_[i];
}

double & psstVariable::time_step()
{
	return time_step_;
}

const double & psstVariable::time_step() const
{
	return time_step_;
}

psstVariable::~psstVariable()
{}

void psstVariable::private_set_dimension(size_t dim)
{
	if (dim > 3)
		throw std::length_error("Error: Invalid dimension!");
	dimension_ = dim;
}
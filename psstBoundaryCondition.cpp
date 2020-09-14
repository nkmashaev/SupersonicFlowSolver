#include "psstVariable.h"
#include "psstBoundaryCondition.h"
#include "psstFiniteVolumeGridTools.h"
#include <vector>
#include <stdexcept>
#include <exception>

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstSupersonicInlet																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstSupersonicInlet::psstSupersonicInlet(const std::vector<double> & velocity, double pressure) :
	velocity_(),
	pressure_(0.0),
	storage_(0)
{
	if (velocity.size() != 2 && velocity.size() != 3)
		throw std::invalid_argument("Error: Wrong dimension!");
	velocity_ = velocity;
	pressure_ = pressure;
}

void psstSupersonicInlet::initialize()
{
	std::shared_ptr<psstInterfaceVariable> var = 0;
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		var = storage_->volumes_.at(i)->variable();

		var->pressure() = pressure_;
		for (size_t j = 0; j < velocity_.size(); ++j)
			var->velocity(j) = velocity_[j];
		calculate_at(i);
	}
}

psstInterfaceBoundaryCondition * psstSupersonicInlet::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstSupersonicInlet(*this);
	return new_bc;
}

void psstSupersonicInlet::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstSupersonicOutlet																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstSupersonicOutlet::initialize()
{
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
		calculate_at(i);
}

psstInterfaceBoundaryCondition * psstSupersonicOutlet::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstSupersonicOutlet(*this);
	return new_bc;
}

void psstSupersonicOutlet::calculate_at(size_t i) const
{
	std::shared_ptr<psstInterfaceVariable> var = storage_->volumes_.at(i)->variable();
	std::shared_ptr<psstInterfaceVariable> neighbour = storage_->volumes_.at(i)->neighbour(0)->variable();

	for (size_t j = 0; j < var->dimension(); ++j)
		var->velocity(j) = neighbour->velocity(j);
	var->pressure() = neighbour->pressure();
}

void psstSupersonicOutlet::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstInviscidWall 																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstInviscidWall::initialize()
{
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
		calculate_at(i);
}

psstInterfaceBoundaryCondition * psstInviscidWall::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstInviscidWall(*this);
	return new_bc;
}

void psstInviscidWall::calculate_at(size_t i) const
{
	double normal_velocity = 0.0;
	std::shared_ptr<psstInterfaceVariable> var = storage_->volumes_.at(i)->variable();
	std::shared_ptr<psstInterfaceVariable> neighbour = storage_->volumes_.at(i)->neighbour(0)->variable();
	std::shared_ptr<const psstInterfaceVector> normal = storage_->volumes_.at(i)->face(0)->normal_vector();
	size_t dim = var->dimension();

	normal_velocity = 0.0;
	for (size_t j = 0; j < dim; ++j)
	{
		normal_velocity += neighbour->velocity(j) * normal->component(j);
	}

	for (size_t j = 0; j < dim; ++j)
	{
		var->velocity(j) = neighbour->velocity(j) - normal_velocity * normal->component(j);
	}
	var->pressure() = neighbour->pressure();
}

void psstInviscidWall::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstInlet																																															|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstInlet::psstInlet(const std::vector<double> & velocity) :
	velocity_(),
	storage_(0)
{
	if (velocity.size() != 2 || velocity.size() != 3)
		throw std::invalid_argument("Error: Wrong dimension!");
	velocity_ = velocity;
}

void psstInlet::initialize()
{
	std::shared_ptr<psstInterfaceVariable> var = 0;
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		var = storage_->volumes_.at(i)->variable();
		for (size_t j = 0; j < velocity_.size(); ++j)
		{
			var->velocity(j) = velocity_[j];
		}
		calculate_at(i);
	}
}


psstInterfaceBoundaryCondition * psstInlet::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstInlet(*this);
	return new_bc;
}

void psstInlet::calculate_at(size_t i) const
{
	storage_->volumes_.at(i)->variable()->pressure() = storage_->volumes_.at(i)->neighbour(0)->variable()->pressure();
}

void psstInlet::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstOutlet																																															|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstOutlet::psstOutlet(double pressure) :
	pressure_(pressure),
	storage_(0)
{}

void psstOutlet::initialize()
{
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		storage_->volumes_.at(i)->variable()->pressure() = pressure_;
		calculate_at(i);
	}
}

psstInterfaceBoundaryCondition * psstOutlet::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstOutlet(*this);
	return new_bc;
}

void psstOutlet::calculate_at(size_t i) const
{
	std::shared_ptr<psstInterfaceVariable> var = storage_->volumes_.at(i)->variable();
	for (size_t j = 0; j < storage_->volumes_.at(i)->variable()->dimension(); ++j)
	{
		var->velocity(j) = storage_->volumes_.at(i)->neighbour(0)->variable()->velocity(j);
	}
}

void psstOutlet::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstThermalInlet																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstThermalInlet::psstThermalInlet(double temperature) :
	temperature_(temperature),
	storage_(0)
{}

void psstThermalInlet::initialize()
{
	std::shared_ptr<psstInterfaceVariable> var = 0;
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		var = storage_->volumes_.at(i)->variable();
		var->temperature() = temperature_;
		calculate_at(i);
	}
}

psstInterfaceBoundaryCondition * psstThermalInlet::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstThermalInlet(*this);
	return new_bc;
}

void psstThermalInlet::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstThermalOutlet																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstThermalOutlet::initialize()
{
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		calculate_at(i);
	}
}

psstInterfaceBoundaryCondition * psstThermalOutlet::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstThermalOutlet(*this);
	return new_bc;
}

void psstThermalOutlet::calculate_at(size_t i) const
{
	storage_->volumes_.at(i)->variable()->temperature() = storage_->volumes_.at(i)->neighbour(0)->variable()->temperature();
}

void psstThermalOutlet::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstThermalWall_T																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstThermalWall_T::psstThermalWall_T(double temperature) :
	temperature_(temperature),
	storage_(0)
{}

void psstThermalWall_T::initialize()
{
	std::shared_ptr<psstInterfaceVariable> var = 0;
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		var = storage_->volumes_.at(i)->variable();
		var->temperature() = temperature_;
		calculate_at(i);
	}
}

psstInterfaceBoundaryCondition * psstThermalWall_T::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstThermalWall_T(*this);
	return new_bc;
}

void psstThermalWall_T::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceBoundaryCondition																																										|
|	psstThermalWall_Q																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstThermalWall_Q::initialize()
{
	for (size_t i = 0; i < storage_->volumes_.size(); ++i)
	{
		calculate_at(i);
	}
}

psstInterfaceBoundaryCondition * psstThermalWall_Q::clone() const
{
	psstInterfaceBoundaryCondition * new_bc = new psstThermalWall_Q(*this);
	return new_bc;
}

void psstThermalWall_Q::calculate_at(size_t i) const
{
	storage_->volumes_.at(i)->variable()->temperature() = storage_->volumes_.at(i)->neighbour(0)->variable()->temperature();
}

void psstThermalWall_Q::get_data(psstVolumeStorage * storage)
{
	if (storage == 0)
		throw std::invalid_argument("Error: Expected boundary data but null_ptr found!");
	storage_ = storage;
}
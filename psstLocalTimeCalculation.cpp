#include "psstLocalTimeCalculation.h"
#include "psstFiniteVolumeGrid.h"
#include "psstFiniteVolumeGridTools.h"
#include "psstVariable.h"
#include <vector>

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceLocalTime																																												|
|	psstUniformTime																																														|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstUniformTime::psstUniformTime(double velocity) :
	ref_velocity_(velocity),
	CFL_(0.5),
	grid_(0)
{}


void psstUniformTime::initialize()
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");

	double dt = 0.0;
	if (grid_->dimension() == 2)
		dt = CFL() * std::sqrt(grid_->min_vol()) / ref_velocity_;
	else
		dt = CFL() * std::pow(grid_->min_vol(), 1.0 / 3.0) / ref_velocity_;

	std::map<size_t, psstVolumeStorage>::iterator iter = grid_->volume_storage()->begin();
	for (iter; iter != grid_->volume_storage()->end(); ++iter)
	{
		if (iter->second.bc_type_id_ == 1)
		{
			for (size_t j = 0; j < iter->second.volumes_.size(); ++j)
			{
				iter->second.volumes_.at(j)->variable()->time_step() = dt;
			}
		}
	}
}

psstInterfaceLocalTime * psstUniformTime::clone() const
{
	psstInterfaceLocalTime * loc_time = new psstUniformTime(*this);
	return loc_time;
}

void psstUniformTime::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceLocalTime																																												|
|	psstLocalAccTime1 																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstLocalAccTime1::initialize()
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");
	dim_ = grid_->dimension();
	rm_ = grid_->rm();
	gamma_ = grid_->gamma();

	std::map<size_t, psstVolumeStorage>::iterator iter = grid_->volume_storage()->begin();
	std::vector<std::shared_ptr<psstInterfaceVolume>> * vols = 0;
	for (iter; iter != grid_->volume_storage()->end(); ++iter)
	{
		if (iter->second.bc_type_id_ != 1)
			continue;

		vols = &iter->second.volumes_;
		for (size_t i = 0; i < vols->size(); ++i)
			calculate_at(vols->at(i));
	}
}

psstInterfaceLocalTime * psstLocalAccTime1::clone() const
{
	psstInterfaceLocalTime * loc_time = new psstLocalAccTime1(*this);
	return loc_time;
}

void psstLocalAccTime1::calculate_at(std::shared_ptr<psstInterfaceVolume> vol) const
{
	if (vol == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");
	if (vol->is_boundary())
		throw std::invalid_argument("Error: Expected inner zone id!");
	var_ = vol->variable();

	double tau_inviscid;
	double volume = vol->volume();
	double face_square = 0.0;
	double normal_velocity = 0.0;
	double normal_neighbour = 0.0;
	double temp = 0.0;
	double c = std::sqrt(gamma_ * rm_ * var_->temperature());
	double c_neighbour = 0.0;
	double denumerator = 0.0;
	size_t face_numb = vol->number_of_faces();
	

	for (size_t j = 0; j < face_numb; ++j)
	{
		normal_velocity = 0.0;
		curr_face_ = vol->face(j);
		normal_ = curr_face_->normal_vector();
		face_square = curr_face_->square();

		neighbour_var_ = curr_face_->neighbour_for(vol)->variable();
		c_neighbour = (std::sqrt(gamma_ * rm_ *neighbour_var_->temperature()) + c) * 0.5;
		for (size_t k = 0; k < dim_; ++k)
		{
			temp = normal_->component(k);
			normal_velocity += var_->velocity(k) * temp;
			normal_neighbour += neighbour_var_->velocity(k) * temp;
		}
		normal_neighbour = (normal_neighbour + normal_velocity) * 0.5;
		denumerator += (std::abs(normal_neighbour) + c_neighbour) * face_square;
	}

	tau_inviscid = CFL() * volume / denumerator;
	var_->time_step() = tau_inviscid;
}

void psstLocalAccTime1::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceLocalTime																																												|
|	psstLocalAccTime2 																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstLocalAccTime2::initialize()
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");
	dim_ = grid_->dimension();
	rm_ = grid_->rm();
	gamma_ = grid_->gamma();
	s_.resize(dim_);
	for (size_t i = 0; i < dim_; ++i)
		s_[i] = 0.0;

	std::map<size_t, psstVolumeStorage>::iterator iter = grid_->volume_storage()->begin();
	std::vector<std::shared_ptr<psstInterfaceVolume>> * vols = 0;
	for (iter; iter != grid_->volume_storage()->end(); ++iter)
	{
		if (iter->second.bc_type_id_ != 1)
			continue;

		vols = &iter->second.volumes_;
		for (size_t i = 0; i < vols->size(); ++i)
			calculate_at(vols->at(i));
	}
}

psstInterfaceLocalTime * psstLocalAccTime2::clone() const
{
	psstInterfaceLocalTime * loc_time = new psstLocalAccTime2(*this);
	return loc_time;
}

void psstLocalAccTime2::calculate_at(std::shared_ptr<psstInterfaceVolume> vol) const
{
	if (vol == 0)
		throw std::invalid_argument("Error: Expected data but null_ptr found!");
	if (vol->is_boundary())
		throw std::invalid_argument("Error: Expected inner zone id!");
	var_ = vol->variable();

	double tau_inviscid;
	double volume = vol->volume();
	double face_square = 0.0;
	double temp = 0.0;
	double c = std::sqrt(gamma_ * rm_ * var_->temperature());
	double denumerator = 0.0;
	size_t face_numb = vol->number_of_faces();


	for (size_t j = 0; j < face_numb; ++j)
	{
		curr_face_ = vol->face(j);
		normal_ = curr_face_->normal_vector();
		face_square = curr_face_->square();

		for (size_t k = 0; k < dim_; ++k)
		{
			s_[k] += face_square * std::abs(normal_->component(k));
		}
	}

	for (size_t j = 0; j < dim_; ++j)
	{
		denumerator += (std::abs(var_->velocity(j)) + c) * s_[j] * 0.5;
		s_[j] = 0.0;
	}

	
	tau_inviscid = CFL() * volume / denumerator;
	var_->time_step() = tau_inviscid;
}

void psstLocalAccTime2::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}
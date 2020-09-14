#include <memory>
#include "psstValueReconstructor.h"
#include "psstFiniteVolumeGridTools.h"
#include "psstScalarLimiters.h"
#include "psstGradCalculation.h"
#include "psstVariable.h"

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceValueReconstructor																																										|
|	psstFirstOrderReconstructor																																											|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstInterfaceValueReconstructor * psstFirstOrderReconstructor::clone() const
{
	psstInterfaceValueReconstructor * rec = new psstFirstOrderReconstructor(*this);
	return rec;
}

void psstFirstOrderReconstructor::attach_limiter(std::shared_ptr<const psstInterfaceScalarLimiter> limiter)
{
	if (limiter.get() == 0)
		throw std::invalid_argument("Error: Expected limiter but null_ptr found!");
	limiter_ = limiter;
}

void psstFirstOrderReconstructor::reconstruct_vals(std::shared_ptr<const psstInterfaceFace> face, std::vector<std::shared_ptr<psstInterfaceVariable>> & rvars)
{
	std::shared_ptr<const psstInterfaceVolume> vols[2];
	std::shared_ptr<const psstInterfaceVariable> vars[2];

	for (size_t j = 0; j < 2; ++j)
	{
		vols[j] = face->neighbour(j);
		vars[j] = vols[j]->variable();

		rvars[j]->pressure() = vars[j]->pressure();
		rvars[j]->temperature() = vars[j]->temperature();
		rvars[j]->density() = vars[j]->density();
		for (size_t i = 0; i < rvars[j]->dimension(); ++i)
		{
			rvars[j]->velocity(i) = vars[j]->velocity(i);
		}
	}


}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceValueReconstructor																																										|
|	psstReconstructorDM																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstInterfaceValueReconstructor * psstReconstructorDM::clone() const
{
	psstInterfaceValueReconstructor * rec = new psstReconstructorDM(*this);
	return rec;
}

void psstReconstructorDM::attach_limiter(std::shared_ptr<const psstInterfaceScalarLimiter> limiter)
{
	if (limiter.get() == 0)
		throw std::invalid_argument("Error: Expected limiter but null_ptr found!");
	limiter_ = limiter;
}

void psstReconstructorDM::reconstruct_vals(std::shared_ptr<const psstInterfaceFace> face, std::vector<std::shared_ptr<psstInterfaceVariable>> & rvars)
{
	std::shared_ptr<const psstInterfaceVolume> vols[2];
	std::shared_ptr<const psstInterfaceVariable> vars[2];
	std::shared_ptr<const psstInterfaceVertex> centers_[2];

	for (size_t j = 0; j < 2; ++j)
	{
		vols[j] = face->neighbour(j);
		vars[j] = vols[j]->variable();
		centers_[j] = vols[j]->center_vertex();

		rvars[j]->pressure() = vars[j]->pressure();
		ptr_vars_[j][dim_] = &rvars[j]->pressure();

		ptr_grads_[j][dim_] = vars[j]->grad_pressure();

		rvars[j]->temperature() = vars[j]->temperature();
		ptr_vars_[j][dim_ + 1] = &rvars[j]->temperature();

		ptr_grads_[j][dim_ + 1] = vars[j]->grad_temperature();

		rvars[j]->density() = vars[j]->density();
		ptr_vars_[j][dim_ + 2] = &rvars[j]->density();

		//ptr_grads_[j][dim_ + 2] = vars[j]->grad_density();
		for (size_t i = 0; i < dim_; ++i)
		{
			rvars[j]->velocity(i) = vars[j]->velocity(i);
			ptr_vars_[j][i] = &rvars[j]->velocity(i);
			ptr_grads_[j][i] = vars[j]->grad_velocity(i);
		}
	}
	
	for (size_t j = 0; j < dim_; ++j)
	{
		R_LR_[j] = centers_[0]->component(j) - centers_[1]->component(j);
	}

	double dot_product_left;
	double dot_product_right;
	for (size_t i = 0; i < dim_ + 2; ++i)
	{
		dot_product_left = 0.0;
		dot_product_right = 0.0;
		for (size_t k = 0; k < dim_; ++k)
		{
			dot_product_left += ptr_grads_[1][i]->component(k) * R_LR_[k];
			dot_product_right += ptr_grads_[0][i]->component(k) * R_LR_[k];
		}
		double val_left_minus = *ptr_vars_[1][i] - (2.0 * dot_product_left - (*ptr_vars_[0][i] - *ptr_vars_[1][i]));
		double val_left_plus = *ptr_vars_[0][i];
		double val_right_minus = *ptr_vars_[1][i];
		double val_right_plus = *ptr_vars_[0][i] + (2.0 * dot_product_right - (*ptr_vars_[0][i] - *ptr_vars_[1][i]));
		double numerator_left = val_left_plus - *ptr_vars_[1][i];
		double denumerator_left = *ptr_vars_[1][i] - val_left_minus;
		double numerator_right = val_right_minus - *ptr_vars_[0][i];
		double denumerator_right = *ptr_vars_[0][i] - val_right_plus;
		*ptr_vars_[0][i] += 0.5 * limiter_->get_value(numerator_right, denumerator_right) * denumerator_right;
		*ptr_vars_[1][i] += 0.5 * limiter_->get_value(numerator_left, denumerator_left) * denumerator_left;
	}

	for (size_t j = 0; j < 2; ++j)
		*ptr_vars_[j][dim_ + 2] = *ptr_vars_[j][dim_] / (rm_ * *ptr_vars_[j][dim_ + 1]);
}
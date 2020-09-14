#include "psstGradCalculation.h"
#include <memory>
#include <map>
#include <vector>
#include "psstFiniteVolumeGridTools.h"
#include "psstVariable.h"
#include "psstFiniteVolumeGrid.h"


/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceGradientCalc																																											|
|	psstGreenGauss																																														|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstGreenGauss::initialize()
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");
	dim_ = grid_->dimension();

	ptr_vars_ = std::vector<std::vector<double *>>(2, std::vector<double*>(dim_ + 3, 0));
	ptr_grads_ = std::vector<std::vector<std::shared_ptr<psstInterfaceVector>>>(2, std::vector<std::shared_ptr<psstInterfaceVector>>(dim_ + 3, 0));
	temp_vals_.resize(dim_ + 3,0.0);

	distances_.clear();
	beta_.clear();
	distances_.resize(2);
	for (size_t i = 0; i < distances_.size(); ++i)
		distances_[i].resize(grid_->number_of_faces());
	beta_.resize(grid_->number_of_faces());
	
	std::shared_ptr<const psstInterfaceVertex> curr_center = 0;
	std::shared_ptr<const psstInterfaceVertex> n_center = 0;

	vol_[0] = 0;
	curr_face_ = 0;
	std::map<size_t, psstFaceStorage>::const_iterator iter = grid_->face_storage()->cbegin();
	double add = 0.0;
	double distance = 0.0;
	for (iter; iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;

		for (size_t i = 0; i < iter->second.faces_.size(); ++i)
		{
			curr_face_ = iter->second.faces_.at(i);
			curr_center = curr_face_->center_vertex();
			for (size_t j = 0; j < 2; ++j)
			{
				distances_[j][i] = 0.0;
				vol_[j] = curr_face_->neighbour(j);				
				if (!vol_[j]->is_boundary())
				{
					n_center = vol_[j]->center_vertex();
					distance = 0.0;
					for (size_t k = 0; k < dim_; ++k)
					{
						add = (curr_center->component(k) - n_center->component(k));
						add = add * add;
						distance += add;
					}
					distances_[j][i] = std::sqrt(distance);
				}
			}

			beta_[i] = distances_[1][i] / (distances_[0][i] + distances_[1][i]);
		}
	}
}

psstInterfaceGradientCalc * psstGreenGauss::clone() const
{
	psstInterfaceGradientCalc * new_grad_calc = new psstGreenGauss(*this);
	return new_grad_calc;
}

void psstGreenGauss::calculate() const
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");

	std::vector<double> velocity_(grid_->dimension(), 0.0);

	double rbeta = 0;
	double face_sqr = 0.0;
	
	std::map<size_t, psstFaceStorage>::const_iterator iter = grid_->face_storage()->cbegin();
	int sign = 1;

	for (iter = grid_->face_storage()->cbegin(); iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;

		for (size_t i = 0; i < iter->second.faces_.size(); ++i)
		{
			curr_face_ = iter->second.faces_.at(i);
			for (size_t j = 0; j < 2; ++j)
			{
				vol_[j] = curr_face_->neighbour(j);
				var_[j] = vol_[j]->variable();
				var_[j]->nullify_grad();
			}
		}
	}

	iter = grid_->face_storage()->cbegin();
	for (iter; iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;

		for (size_t i = 0; i < iter->second.faces_.size(); ++i)
		{
			curr_face_ = iter->second.faces_.at(i);
			normal_ = curr_face_->normal_vector();
			face_sqr = curr_face_->square();
			for (size_t j = 0; j < 2; ++j)
			{
				vol_[j] = curr_face_->neighbour(j);
				var_[j] = vol_[j]->variable();
				geom_par_[j] = face_sqr / vol_[j]->volume();
				ptr_vars_[j][dim_] = &var_[j]->pressure();
				ptr_grads_[j][dim_] = var_[j]->grad_pressure();
				ptr_vars_[j][dim_ + 1] = &var_[j]->temperature();
				ptr_grads_[j][dim_ + 1] = var_[j]->grad_temperature();
				ptr_vars_[j][dim_ + 2] = &var_[j]->density();
				ptr_grads_[j][dim_ + 2] = var_[j]->grad_density();
				for (size_t n = 0; n < dim_; ++n)
				{
					ptr_vars_[j][n] = &var_[j]->velocity(n);
					ptr_grads_[j][n] = var_[j]->grad_velocity(n);
				}
			}

			rbeta = 1 - beta_[i];
			for (size_t n = 0; n < temp_vals_.size(); ++n)
			{
				temp_vals_[n] = beta_[i] * *ptr_vars_[0][n] + rbeta * *ptr_vars_[1][n];
				sign = 1;
				for (size_t j = 0; j < 2; ++j)
				{
					sign *= -1;
					for (size_t k = 0; k < dim_; ++k)
						ptr_grads_[j][n]->component(k) += sign * temp_vals_[n] * geom_par_[j] * normal_->component(k);
				}
			}
		}
	}

	for (iter = grid_->face_storage()->cbegin(); iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;

		for (size_t i = 0; i < iter->second.faces_.size(); ++i)
		{
			curr_face_ = iter->second.faces_.at(i);
			for (size_t j = 0; j < 2; ++j)
			{
				vol_[j] = curr_face_->neighbour(j);
				var_[j] = vol_[j]->variable();
				var_[j]->grad_pressure()->initialize();
				var_[j]->grad_density()->initialize();
				var_[j]->grad_temperature()->initialize();
				for (size_t n = 0; n < dim_; ++n)
				{
					var_[j]->grad_velocity(n)->initialize();
				}
			}
		}
	}

}

void psstGreenGauss::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceGradientCalc																																											|
|	psstLeastSquare																																														|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstLeastSquare::initialize()
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");
	

	size_t dim = grid_->dimension();
	size_t vol_numb = grid_->number_of_inner_volumes();
	size_t coeff_size = (dim - 1) * dim / 2 + dim;
	std::vector<double> coeffs(coeff_size, 0.0);
	beta_.resize(vol_numb);
	r_ = std::vector<std::vector<double>>(vol_numb, std::vector<double>(coeff_size, 0.0));
	ptr_grads_.resize(dim + 3);
	std::vector<double> diff(dim, 0.0);
	double w = 0.0;
	std::map<size_t, psstVolumeStorage>::iterator iter = grid_->volume_storage()->begin();
	std::shared_ptr<psstInterfaceVolume> vol = 0;
	std::shared_ptr<psstInterfaceVolume> neighbour = 0;
	std::shared_ptr<const psstInterfaceVertex> vol_center = 0;
	std::shared_ptr<const psstInterfaceVertex> neighbour_center = 0;

	double w_sqr;
	for (iter; iter != grid_->volume_storage()->end(); ++iter)
	{
		if (iter->second.bc_type_id_ != 1)
			continue;

		for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
		{
			vol = iter->second.volumes_.at(i);
			vol_center = vol->center_vertex();
			for (size_t k = 0; k < coeffs.size(); ++k)
				coeffs[k] = 0.0;

			for (size_t j = 0; j < vol->number_of_faces(); ++j)
			{
				neighbour = vol->neighbour(j);
				neighbour_center = neighbour->center_vertex();
				w = 0.0;
				for (size_t k = 0; k < dim; ++k)
				{
					diff[k] = neighbour_center->component(k) - vol_center->component(k);
					w += diff[k] * diff[k];
				}
				w_sqr = 1.0 / w;
				coeffs[a] += w_sqr * diff[0] * diff[0];
				coeffs[b] += w_sqr * diff[0] * diff[1];
				coeffs[d] += w_sqr * diff[1] * diff[1];
				if (dim == 3)
				{
					coeffs[c] += w_sqr * diff[0] * diff[2];
					coeffs[e] += w_sqr * diff[1] * diff[2];
					coeffs[g] += w_sqr * diff[2] * diff[2];
				}
			}

			r_[i][a] = std::sqrt(coeffs[a]);
			r_[i][b] = coeffs[b]/r_[i][a];
			r_[i][d] = std::sqrt(coeffs[d] - r_[i][b] * r_[i][b]);
			beta_[i] = 0.0;
			if (dim == 3)
			{
				r_[i][c] = coeffs[c] / r_[i][a];
				r_[i][e] = (coeffs[e] - r_[i][b]/r_[i][a] * coeffs[c]) / r_[i][d];
				r_[i][g] = std::sqrt(coeffs[g] - (r_[i][c] * r_[i][c] + r_[i][e] * r_[i][e]));
				beta_[i] = (r_[i][b] * r_[i][e] - r_[i][c] * r_[i][d]) / (r_[i][a] * r_[i][d]);
			}
		}
	}
}

psstInterfaceGradientCalc * psstLeastSquare::clone() const
{
	psstInterfaceGradientCalc * new_grad_calc = new psstLeastSquare(*this);
	return new_grad_calc;
}

void psstLeastSquare::calculate() const
{
	if (grid_ == 0)
		throw std::exception("Error: Grid data is not attached!");

	std::map<size_t, psstVolumeStorage>::iterator iter = grid_->volume_storage()->begin();
	std::shared_ptr<psstInterfaceVolume> vol = 0;
	std::shared_ptr<psstInterfaceVolume> neighbour = 0;
	std::shared_ptr<const psstInterfaceVertex> vol_center = 0;
	std::shared_ptr<const psstInterfaceVertex> neighbour_center = 0;
	std::shared_ptr<psstInterfaceVariable> var = 0;
	std::shared_ptr<psstInterfaceVariable> neighbour_var = 0;
	size_t dim = grid_->dimension();
	std::vector<double> diff(dim, 0.0);
	double w = 0.0;
	double w_sqr = 0.0;
	std::vector<double> alpha(dim, 0.0);
	std::vector<double> theta(dim, 0.0);
	std::vector<std::vector<double>> phi(dim + 3, std::vector<double>(dim, 0.0));

	for (iter; iter != grid_->volume_storage()->end(); ++iter)
	{
		if (iter->second.bc_type_id_ != 1)
			continue;
		
		for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
		{
			vol = iter->second.volumes_.at(i);
			vol_center = vol->center_vertex();
			var = vol->variable();
			var->nullify_grad();

			for (size_t k = 0; k < dim; ++k)
			{
				ptr_grads_[k] = var->grad_velocity(k);
			}
			ptr_grads_[dim] = var->grad_pressure();
			ptr_grads_[dim + 1] = var->grad_temperature();
			ptr_grads_[dim + 2] = var->grad_density();


			for (size_t j = 0; j < vol->number_of_faces(); ++j)
			{
				neighbour = vol->neighbour(j);
				neighbour_var = neighbour->variable();
				neighbour_center = neighbour->center_vertex();
				w = 0.0;
				for (size_t k = 0; k < dim; ++k)
				{
					diff[k] = neighbour_center->component(k) - vol_center->component(k);
					w += diff[k] * diff[k];
					alpha[k] = 0.0;
					for (size_t n = 0; n < dim; ++n)
						phi[n][k] = (neighbour_var->velocity(n) - var->velocity(n));
					phi[dim][k] =  (neighbour_var->pressure() - var->pressure());
					phi[dim + 1][k] = (neighbour_var->temperature() - var->temperature());
					phi[dim + 2][k] = (neighbour_var->density() - var->density());

					
				}


				w_sqr = 1.0 / w;
				if (dim == 2)
				{
					alpha[0] = diff[0] / (r_[i][a] * r_[i][a]);
					alpha[1] = (diff[1] - r_[i][b] / r_[i][a] * diff[0]) / (r_[i][d] * r_[i][d]);
				}
				else
				{
					alpha[3] = (diff[2] - r_[i][e] / r_[i][d] * diff[1] + beta_[i] * diff[2]) / (r_[i][g] * r_[i][g]);
				}

				theta[0] = alpha[0] - r_[i][b] / r_[i][a] * alpha[1];
				theta[1] = alpha[1];
				if (dim == 3)
				{
					theta[0] += alpha[2] * beta_[i];
					theta[1] -= r_[i][e] / r_[i][d] * alpha[2];
					theta[2] = alpha[2];
				}				

				for (size_t n = 0; n < dim + 3; ++n)
				{
					for (size_t k = 0; k < dim; ++k)
					{
						ptr_grads_[n]->component(k) += w_sqr * theta[k] * phi[n][k];
					}
				}
			}
			
			for (size_t n = 0; n < dim + 3; ++n)
			{
				ptr_grads_[n]->initialize();
			}
		}

		
	}

}

void psstLeastSquare::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}
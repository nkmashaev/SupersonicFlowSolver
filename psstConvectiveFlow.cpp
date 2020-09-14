#include "psstConvectiveFlow.h"
#include "psstFiniteVolumeGrid.h"
#include "psstFiniteVolumeGridTools.h"
#include "psstValueReconstructor.h"
#include "psstVariable.h"
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>

double average_by_Roe(double left_val, double right_val, double left_density, double right_density)
{
	double value_averaged_by_Roe = 0.0;
	double numerator = 0.0;
	double denumerator = 0.0;
	double left_d = std::sqrt(left_density);
	double right_d = std::sqrt(right_density);
	numerator = left_val * left_d + right_val * right_d;
	denumerator = (left_d + right_d);
	value_averaged_by_Roe = numerator / denumerator;
	return value_averaged_by_Roe;
}

double average_by_Roe(double left_density, double right_density)
{
	double densityAveragedByRoe = 0.0;
	densityAveragedByRoe = std::sqrt(left_density * right_density);
	return densityAveragedByRoe;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceConvectiveFlow																																											|
|	psstConvectiveHLL																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstConvectiveHLL::initialize()
{
	face_sqr_ = 0.0;
	dim_ = grid_->dimension();
	vel_.resize(2);
	flow_i_.resize(2);
	w_i_.resize(2);
	v_Roe_.resize(dim_);
	flow_.resize(dim_ + 2);
	rvars_.resize(2);
	a_Roe_ = 0.0;
	H_Roe_ = 0.0;
	v_Roe_n_ = 0.0;
	v_sqr_ = 0.0;
	for (size_t i = 0; i < flow_.size(); ++i)
		flow_[i] = 0.0;

	w_i_.resize(2);
	flow_i_.resize(2);
	for (size_t i = 0; i < 2; ++i)
	{
		pressure_[i] = 0.0;
		temperature_[i] = 0.0;
		density_[i] = 0.0;
		v_n_[i] = 0.0;
		H_[i] = 0.0;
		E_[i] = 0.0;
		a_[i] = 0.0;
		s_[i] = 0.0;
		vols_[i] = 0;
		vars_[i] = 0;
		vel_[i].resize(dim_);
		for (size_t j = 0; j < dim_; ++j)
			vel_[i][j] = 0.0;
		w_i_[i].resize(flow_.size());
		flow_i_[i].resize(flow_.size());
		for (size_t j = 0; j < flow_.size(); ++j)
		{
			w_i_[i][j] = 0.0;
			flow_i_[i][j] = 0.0;
		}
	}
	
	std::map<size_t, psstVolumeStorage>::const_iterator iter = grid_->volume_storage()->cbegin();

	psstInterfaceVariable * var = iter->second.volumes_[0]->variable()->clone();
	rvars_[0].reset(var);

	var = iter->second.volumes_[0]->variable()->clone();
	rvars_[1].reset(var);

	curr_face_ = 0;
	curr_normal_ = 0;
}

psstInterfaceConvectiveFlow * psstConvectiveHLL::clone() const
{
	psstInterfaceConvectiveFlow * convective = new psstConvectiveHLL(*this);
	return convective;
}

void psstConvectiveHLL::attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor> reconstructor)
{
	reconstructor_ = reconstructor;
}

void psstConvectiveHLL::calculate()
{
	std::map<size_t, psstFaceStorage>::const_iterator iter = grid_->face_storage()->cbegin();
	size_t dim_plus = dim_ + 1;
	size_t w_size = dim_ + 2;
	double cp = grid_->cp();
	double cv = grid_->cv();
	int sign = 1;
	double rm = grid_->rm();
	double gamma = grid_->gamma();
	double temp = 0.0;

	for (iter; iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;
 
		if (iter->second.bc_type_id_ != 2)
		{
			
			for (size_t i = 0; i < iter->second.faces_.size(); ++i)
			{
			
				/*out_file << iter->first << " " << i << " ";*/
				curr_face_ = iter->second.faces_.at(i);
				face_sqr_ = curr_face_->square();
				curr_normal_ = curr_face_->normal_vector();
				for (size_t j = 0; j < 2; ++j)
				{
					vols_[j] = curr_face_->neighbour(j);
					vars_[j] = vols_[j]->variable();
				}
				
				if (vols_[0]->is_boundary())
				{

					v_n_[0] = 0.0;
					v_sqr_ = 0.0;
					pressure_[0] = vars_[0]->pressure();
					temperature_[0] = vars_[0]->temperature();
					density_[0] = vars_[0]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[0][k] = vars_[0]->velocity(k);
						v_n_[0] += vel_[0][k] * curr_normal_->component(k);
						v_sqr_ += vel_[0][k] * vel_[0][k];
					}

					temperature_[0] = vars_[0]->temperature();
					H_[0] = cp * temperature_[0] + 0.5 * v_sqr_;

					for (size_t k = 0; k < dim_; ++k)
					{
						flow_[k] = density_[0] * vel_[0][k] * v_n_[0]
							+ pressure_[0] * curr_normal_->component(k);
					}
					flow_[dim_] = density_[0] * v_n_[0];
					flow_[dim_plus] = density_[0] * v_n_[0] * H_[0];
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[1]->residual(k) += flow_[k] * face_sqr_;
						//out_file << 1 << " " << flow[k] << " " << face_sqr << " " << vars[1]->residual(k) << "\n";
					}
				}
				else 
				{
					v_n_[1] = 0.0;
					v_sqr_ = 0.0;
					pressure_[1] = vars_[1]->pressure();
					temperature_[1] = vars_[1]->temperature();
					density_[1] = vars_[1]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[1][k] = vars_[1]->velocity(k);
						v_n_[1] += vel_[1][k] * curr_normal_->component(k);
						v_sqr_ += vel_[1][k] * vel_[1][k];
					}

					temperature_[1] = vars_[1]->temperature();
					H_[1] = cp * temperature_[1] + 0.5 * v_sqr_;

					for (size_t k = 0; k < dim_; ++k)
					{
						flow_[k] = density_[1] * vel_[1][k] * v_n_[1]
							+ pressure_[1] * curr_normal_->component(k);
					}
					flow_[dim_] = density_[1] * v_n_[1];
					flow_[dim_plus] = density_[1] * v_n_[1] * H_[1];
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[0]->residual(k) -= flow_[k] * face_sqr_;
						//out_file << 1 << " " << flow[k] << " " << face_sqr << " " << vars[1]->residual(k) << "\n";
					}
				}
		
			}
		}
		else
		{
			for (size_t i = 0; i < iter->second.faces_.size(); ++i)
			{
				/*out_file << i << std::endl;*/
				curr_face_ = iter->second.faces_.at(i);
				face_sqr_ = curr_face_->square();
				curr_normal_ = curr_face_->normal_vector();

				reconstructor_->reconstruct_vals(curr_face_, rvars_);
				for (size_t j = 0; j < 2; ++j)
				{
					vols_[j] = curr_face_->neighbour(j);
					vars_[j] = vols_[j]->variable();
					v_n_[j] = 0.0;
					v_sqr_ = 0.0;
					pressure_[j] = rvars_[j]->pressure();
					temperature_[j] = rvars_[j]->temperature();
					density_[j] = rvars_[j]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[j][k] = rvars_[j]->velocity(k);
						v_n_[j] += vel_[j][k] * curr_normal_->component(k);
						v_sqr_ += vel_[j][k] * vel_[j][k];
					}

					a_[j] = std::sqrt(gamma * rm * temperature_[j]);
					H_[j] = cp * temperature_[j] + 0.5 * v_sqr_;
					E_[j] = cv * temperature_[j] + 0.5 * v_sqr_;
					for (size_t k = 0; k < dim_; ++k)
					{
						flow_i_[j][k] = density_[j] * vel_[j][k] * v_n_[j]
							+ pressure_[j] * curr_normal_->component(k);
						w_i_[j][k] = density_[j] * vel_[j][k];
					}
					w_i_[j][dim_] = density_[j];
					w_i_[j][dim_plus] = density_[j] * E_[j];
					flow_i_[j][dim_] = density_[j] * v_n_[j];
					flow_i_[j][dim_plus] = density_[j] * v_n_[j] * H_[j];
				}

				H_Roe_ = average_by_Roe(H_[1], H_[0], density_[1], density_[0]);
				v_Roe_n_= 0.0;
				v_sqr_ = 0.0;
				for (size_t k = 0; k < dim_; ++k)
				{
					v_Roe_[k] = average_by_Roe(vel_[1][k], vel_[0][k], density_[1], density_[0]);
					v_Roe_n_ += v_Roe_[k] * curr_normal_->component(k);
					v_sqr_ += v_Roe_[k] * v_Roe_[k];
				}
				a_Roe_ = std::sqrt((gamma - 1.0) * (H_Roe_ - 0.5 * v_sqr_));

				s_[1] = v_n_[1] - a_[1];
				temp = v_Roe_n_ - a_Roe_;
				if (s_[1] > temp)
					s_[1] = temp;

				s_[0] = v_n_[0] + a_[0];
				temp = v_Roe_n_ + a_Roe_;
				if (s_[0] < temp)
					s_[0] = temp;

				if (s_[1] >= 0.0)
				{
					for (size_t j = 0; j < w_size; ++j)
						flow_[j] = flow_i_[1][j];
				}
				else if (s_[1] <= 0.0 && s_[0] >= 0.0)
				{
					for (size_t j = 0; j < w_size; ++j)
						flow_[j] = (s_[0] * flow_i_[1][j] - s_[1] * flow_i_[0][j] + s_[0] * s_[1] * (w_i_[0][j] - w_i_[1][j])) / (s_[0] - s_[1]);
				}
				else if (s_[0] <= 0.0)
				{
					for (size_t j = 0; j < w_size; ++j)
						flow_[j] = flow_i_[0][j];
				}
				else
				{
					throw std::exception("Error: Could not calculate convective HLL flow!");
				}

				/*out_file << iter->first << " " << i << " " << H_Roe << " " << a_Roe << " " << v_Roe_normal << " "<< v_Roe_sqr << " ";*/
				sign = 1;
				for (size_t j = 0; j < 2; ++j)
				{
					sign = (-1) * sign;
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[j]->residual(k) += sign * flow_[k] * face_sqr_;
						/*out_file << sign << " " << flow[k] << " " << face_sqr << " " << vars[j]->residual(k) << std::endl;*/
					}
					//out_file << sign << " " << H[j] << a[j] << normal_velocity[j] << " " << v_sqr[j] << " ";
				}

				//for (size_t j = 0; j < w_size; ++j)
				//	out_file << flow[j] << " ";
				//out_file << std::endl;
			}
		}
	}

	/*out_file.close();*/
}

void psstConvectiveHLL::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceConvectiveFlow																																											|
|	psstConvectiveHLLC																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstConvectiveHLLC::initialize()
{
	face_sqr_ = 0.0;
	dim_ = grid_->dimension();
	vel_.resize(2);
	flow_i_.resize(2);
	w_i_.resize(2);
	v_Roe_.resize(dim_);
	flow_.resize(dim_ + 2);
	rvars_.resize(2);
	a_Roe_ = 0.0;
	H_Roe_ = 0.0;
	v_Roe_n_ = 0.0;
	v_sqr_ = 0.0;
	for (size_t i = 0; i < flow_.size(); ++i)
		flow_[i] = 0.0;

	s_hllc_ = 0.0;
	w_i_.resize(2);
	flow_i_.resize(2);
	for (size_t i = 0; i < 2; ++i)
	{
		pressure_[i] = 0.0;
		temperature_[i] = 0.0;
		density_[i] = 0.0;
		v_n_[i] = 0.0;
		H_[i] = 0.0;
		E_[i] = 0.0;
		a_[i] = 0.0;
		s_[i] = 0.0;
		vols_[i] = 0;
		vars_[i] = 0;
		vel_[i].resize(dim_);
		for (size_t j = 0; j < dim_; ++j)
			vel_[i][j] = 0.0;
		w_i_[i].resize(flow_.size());
		flow_i_[i].resize(flow_.size());
		for (size_t j = 0; j < flow_.size(); ++j)
		{
			w_i_[i][j] = 0.0;
			flow_i_[i][j] = 0.0;
		}
	}

	std::map<size_t, psstVolumeStorage>::const_iterator iter = grid_->volume_storage()->cbegin();

	psstInterfaceVariable * var = iter->second.volumes_[0]->variable()->clone();
	rvars_[0].reset(var);

	var = iter->second.volumes_[0]->variable()->clone();
	rvars_[1].reset(var);

	curr_face_ = 0;
	curr_normal_ = 0;
}

psstInterfaceConvectiveFlow * psstConvectiveHLLC::clone() const
{
	psstInterfaceConvectiveFlow * convective = new psstConvectiveHLLC(*this);
	return convective;
}

void psstConvectiveHLLC::attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor> reconstructor)
{
	reconstructor_ = reconstructor;
}

void psstConvectiveHLLC::calculate() 
{
	std::map<size_t, psstFaceStorage>::const_iterator iter = grid_->face_storage()->cbegin();
	size_t dim_plus = dim_ + 1;
	size_t w_size = dim_ + 2;
	double cp = grid_->cp();
	double cv = grid_->cv();
	int sign = 1;
	double rm = grid_->rm();
	double gamma = grid_->gamma();
	double temp = 0.0;

	for (iter; iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;

		if (iter->second.bc_type_id_ != 2)
		{

			for (size_t i = 0; i < iter->second.faces_.size(); ++i)
			{

				/*out_file << iter->first << " " << i << " ";*/
				curr_face_ = iter->second.faces_.at(i);
				face_sqr_ = curr_face_->square();
				curr_normal_ = curr_face_->normal_vector();
				for (size_t j = 0; j < 2; ++j)
				{
					vols_[j] = curr_face_->neighbour(j);
					vars_[j] = vols_[j]->variable();
				}

				if (vols_[0]->is_boundary())
				{

					v_n_[0] = 0.0;
					v_sqr_ = 0.0;
					pressure_[0] = vars_[0]->pressure();
					temperature_[0] = vars_[0]->temperature();
					density_[0] = vars_[0]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[0][k] = vars_[0]->velocity(k);
						v_n_[0] += vel_[0][k] * curr_normal_->component(k);
						v_sqr_ += vel_[0][k] * vel_[0][k];
					}

					temperature_[0] = vars_[0]->temperature();
					H_[0] = cp * temperature_[0] + 0.5 * v_sqr_;

					for (size_t k = 0; k < dim_; ++k)
					{
						flow_[k] = density_[0] * vel_[0][k] * v_n_[0]
							+ pressure_[0] * curr_normal_->component(k);
					}
					flow_[dim_] = density_[0] * v_n_[0];
					flow_[dim_plus] = density_[0] * v_n_[0] * H_[0];
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[1]->residual(k) += flow_[k] * face_sqr_;
						//out_file << 1 << " " << flow[k] << " " << face_sqr << " " << vars[1]->residual(k) << "\n";
					}
				}
				else
				{
					v_n_[1] = 0.0;
					v_sqr_ = 0.0;
					pressure_[1] = vars_[1]->pressure();
					temperature_[1] = vars_[1]->temperature();
					density_[1] = vars_[1]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[1][k] = vars_[1]->velocity(k);
						v_n_[1] += vel_[1][k] * curr_normal_->component(k);
						v_sqr_ += vel_[1][k] * vel_[1][k];
					}

					temperature_[1] = vars_[1]->temperature();
					H_[1] = cp * temperature_[1] + 0.5 * v_sqr_;

					for (size_t k = 0; k < dim_; ++k)
					{
						flow_[k] = density_[1] * vel_[1][k] * v_n_[1]
							+ pressure_[1] * curr_normal_->component(k);
					}
					flow_[dim_] = density_[1] * v_n_[1];
					flow_[dim_plus] = density_[1] * v_n_[1] * H_[1];
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[0]->residual(k) -= flow_[k] * face_sqr_;
						//out_file << 1 << " " << flow[k] << " " << face_sqr << " " << vars[1]->residual(k) << "\n";
					}
				}

			}
		}
		else
		{
			for (size_t i = 0; i < iter->second.faces_.size(); ++i)
			{
				/*out_file << i << std::endl;*/
				curr_face_ = iter->second.faces_.at(i);
				face_sqr_ = curr_face_->square();
				curr_normal_ = curr_face_->normal_vector();

				reconstructor_->reconstruct_vals(curr_face_, rvars_);
				for (size_t j = 0; j < 2; ++j)
				{
					vols_[j] = curr_face_->neighbour(j);
					vars_[j] = vols_[j]->variable();
					v_n_[j] = 0.0;
					v_sqr_ = 0.0;
					pressure_[j] = rvars_[j]->pressure();
					temperature_[j] = rvars_[j]->temperature();
					density_[j] = rvars_[j]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[j][k] = rvars_[j]->velocity(k);
						v_n_[j] += vel_[j][k] * curr_normal_->component(k);
						v_sqr_ += vel_[j][k] * vel_[j][k];
					}

					a_[j] = std::sqrt(gamma * rm * temperature_[j]);
					H_[j] = cp * temperature_[j] + 0.5 * v_sqr_;
					E_[j] = cv * temperature_[j] + 0.5 * v_sqr_;
					for (size_t k = 0; k < dim_; ++k)
					{
						flow_i_[j][k] = density_[j] * vel_[j][k] * v_n_[j]
							+ pressure_[j] * curr_normal_->component(k);
						w_i_[j][k] = density_[j] * vel_[j][k];
					}
					w_i_[j][dim_] = density_[j];
					w_i_[j][dim_plus] = density_[j] * E_[j];
					flow_i_[j][dim_] = density_[j] * v_n_[j];
					flow_i_[j][dim_plus] = density_[j] * v_n_[j] * H_[j];
				}

				H_Roe_ = average_by_Roe(H_[1], H_[0], density_[1], density_[0]);
				v_Roe_n_ = 0.0;
				v_sqr_ = 0.0;
				for (size_t k = 0; k < dim_; ++k)
				{
					v_Roe_[k] = average_by_Roe(vel_[1][k], vel_[0][k], density_[1], density_[0]);
					v_Roe_n_ += v_Roe_[k] * curr_normal_->component(k);
					v_sqr_ += v_Roe_[k] * v_Roe_[k];
				}
				a_Roe_ = std::sqrt((gamma - 1.0) * (H_Roe_ - 0.5 * v_sqr_));

				s_[1] = v_n_[1] - a_[1];
				temp = v_Roe_n_ - a_Roe_;
				if (s_[1] > temp)
					s_[1] = temp;

				s_[0] = v_n_[0] + a_[0];
				temp = v_Roe_n_ + a_Roe_;
				if (s_[0] < temp)
					s_[0] = temp;

				s_hllc_ = (pressure_[0] - pressure_[1] + density_[1] * v_n_[1] * (s_[1] - v_n_[1]) - density_[0] * v_n_[0] * (s_[0] - v_n_[0])) / (density_[1] * (s_[1] - v_n_[1]) - density_[0] * (s_[0] - v_n_[0]));

				if (s_[1] >= 0.0)
				{
					for (size_t j = 0; j < w_size; ++j)
						flow_[j] = flow_i_[1][j];
				}
				else if (s_[1] <= 0.0 && s_hllc_ >= 0.0)
				{
					sign = 1;
					temp = s_[sign] - s_hllc_;
					for (size_t j = 0; j < dim_; ++j)
						flow_[j] = (s_hllc_ * (s_[sign] * w_i_[sign][j] - flow_i_[sign][j]) + s_[sign] * (pressure_[sign] + density_[sign] * (s_[sign] - v_n_[sign]) * (s_hllc_ - v_n_[sign])) * curr_normal_->component(j)) / temp;
					flow_[dim_] = (s_hllc_ * (s_[sign] * w_i_[sign][dim_] - flow_i_[sign][dim_])) / temp;
					flow_[dim_plus] = (s_hllc_ * (s_[sign] * w_i_[sign][dim_plus] - flow_i_[sign][dim_plus]) + s_[sign] * (pressure_[sign] + density_[sign] * (s_[sign] - v_n_[sign]) * (s_hllc_ - v_n_[sign])) * s_hllc_) / temp;
				}
				else if (s_hllc_ <= 0.0 && s_[0] >= 0.0)
				{
					sign = 0;
					temp = s_[sign] - s_hllc_;
					for (size_t j = 0; j < dim_; ++j)
						flow_[j] = (s_hllc_ * (s_[sign] * w_i_[sign][j] - flow_i_[sign][j]) + s_[sign] * (pressure_[sign] + density_[sign] * (s_[sign] - v_n_[sign]) * (s_hllc_ - v_n_[sign])) * curr_normal_->component(j)) / temp;
					flow_[dim_] = (s_hllc_ * (s_[sign] * w_i_[sign][dim_] - flow_i_[sign][dim_])) / temp;
					flow_[dim_plus] = (s_hllc_ * (s_[sign] * w_i_[sign][dim_plus] - flow_i_[sign][dim_plus]) + s_[sign] * (pressure_[sign] + density_[sign] * (s_[sign] - v_n_[sign]) * (s_hllc_ - v_n_[sign])) * s_hllc_) / temp;
				}
				else if (s_[0] <= 0.0)
				{
					for (size_t j = 0; j < w_size; ++j)
						flow_[j] = flow_i_[0][j];
				}
				else
				{
					throw std::exception("Error: Could not calculate convective HLLC flow!");
				}

				/*out_file << iter->first << " " << i << " " << H_Roe << " " << a_Roe << " " << v_Roe_normal << " "<< v_Roe_sqr << " ";*/
				sign = 1;
				for (size_t j = 0; j < 2; ++j)
				{
					sign = (-1) * sign;
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[j]->residual(k) += sign * flow_[k] * face_sqr_;
						/*out_file << sign << " " << flow[k] << " " << face_sqr << " " << vars[j]->residual(k) << std::endl;*/
					}
					//out_file << sign << " " << H[j] << a[j] << normal_velocity[j] << " " << v_sqr[j] << " ";
				}

				//for (size_t j = 0; j < w_size; ++j)
				//	out_file << flow[j] << " ";
				//out_file << std::endl;
			}
		}
	}
}

void psstConvectiveHLLC::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceConvectiveFlow																																											|
|	psstConvectiveAUSM																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstConvectiveAUSM::initialize()
{
	face_sqr_ = 0.0;
	dim_ = grid_->dimension();
	vel_.resize(2);
	p_flow_.resize(dim_);
	flow_.resize(dim_ + 2);
	rvars_.resize(2);

	for (size_t i = 0; i < flow_.size(); ++i)
		flow_[i] = 0.0;

	for (size_t i = 0; i < 2; ++i)
	{
		pressure_[i] = 0.0;
		temperature_[i] = 0.0;
		density_[i] = 0.0;
		v_n_[i] = 0.0;
		H_[i] = 0.0;
		M_[i] = 0.0;
		a_[i] = 0.0;
		M_sign_[i] = 0.0;
		p_sign_[i] = 0.0;
		vols_[i] = 0;
		vars_[i] = 0;
		vel_[i].resize(dim_);
		for (size_t j = 0; j < dim_; ++j)
			vel_[i][j] = 0.0;
	}

	std::map<size_t, psstVolumeStorage>::const_iterator iter = grid_->volume_storage()->cbegin();

	psstInterfaceVariable * var = iter->second.volumes_[0]->variable()->clone();
	rvars_[0].reset(var);

	var = iter->second.volumes_[0]->variable()->clone();
	rvars_[1].reset(var);

	curr_face_ = 0;
	curr_normal_ = 0;
}

psstInterfaceConvectiveFlow * psstConvectiveAUSM::clone() const
{
	psstInterfaceConvectiveFlow * convective = new psstConvectiveAUSM(*this);
	return convective;
}

void psstConvectiveAUSM::attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor> reconstructor)
{
	reconstructor_ = reconstructor;
}

void psstConvectiveAUSM::calculate()
{
	std::map<size_t, psstFaceStorage>::const_iterator iter = grid_->face_storage()->cbegin();
	size_t dim_plus = dim_ + 1;
	size_t w_size = dim_ + 2;
	double cp = grid_->cp();
	double cv = grid_->cv();
	int sign = 1;
	double rm = grid_->rm();
	double gamma = grid_->gamma();
	double temp = 0.0;

	for (iter; iter != grid_->face_storage()->cend(); ++iter)
	{
		if (iter->second.bc_type_id_ == 31)
			continue;

		if (iter->second.bc_type_id_ != 2)
		{

			for (size_t i = 0; i < iter->second.faces_.size(); ++i)
			{

				/*out_file << iter->first << " " << i << " ";*/
				curr_face_ = iter->second.faces_.at(i);
				face_sqr_ = curr_face_->square();
				curr_normal_ = curr_face_->normal_vector();
				for (size_t j = 0; j < 2; ++j)
				{
					vols_[j] = curr_face_->neighbour(j);
					vars_[j] = vols_[j]->variable();
				}

				if (vols_[0]->is_boundary())
				{

					v_n_[0] = 0.0;
					v_sqr_ = 0.0;
					pressure_[0] = vars_[0]->pressure();
					temperature_[0] = vars_[0]->temperature();
					density_[0] = vars_[0]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[0][k] = vars_[0]->velocity(k);
						v_n_[0] += vel_[0][k] * curr_normal_->component(k);
						v_sqr_ += vel_[0][k] * vel_[0][k];
					}

					temperature_[0] = vars_[0]->temperature();
					H_[0] = cp * temperature_[0] + 0.5 * v_sqr_;

					for (size_t k = 0; k < dim_; ++k)
					{
						flow_[k] = density_[0] * vel_[0][k] * v_n_[0]
							+ pressure_[0] * curr_normal_->component(k);
					}
					flow_[dim_] = density_[0] * v_n_[0];
					flow_[dim_plus] = density_[0] * v_n_[0] * H_[0];
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[1]->residual(k) += flow_[k] * face_sqr_;
						//out_file << 1 << " " << flow[k] << " " << face_sqr << " " << vars[1]->residual(k) << "\n";
					}
				}
				else
				{
					v_n_[1] = 0.0;
					v_sqr_ = 0.0;
					pressure_[1] = vars_[1]->pressure();
					temperature_[1] = vars_[1]->temperature();
					density_[1] = vars_[1]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[1][k] = vars_[1]->velocity(k);
						v_n_[1] += vel_[1][k] * curr_normal_->component(k);
						v_sqr_ += vel_[1][k] * vel_[1][k];
					}

					temperature_[1] = vars_[1]->temperature();
					H_[1] = cp * temperature_[1] + 0.5 * v_sqr_;

					for (size_t k = 0; k < dim_; ++k)
					{
						flow_[k] = density_[1] * vel_[1][k] * v_n_[1]
							+ pressure_[1] * curr_normal_->component(k);
					}
					flow_[dim_] = density_[1] * v_n_[1];
					flow_[dim_plus] = density_[1] * v_n_[1] * H_[1];
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[0]->residual(k) -= flow_[k] * face_sqr_;
						//out_file << 1 << " " << flow[k] << " " << face_sqr << " " << vars[1]->residual(k) << "\n";
					}
				}

			}
		}
		else
		{
			for (size_t i = 0; i < iter->second.faces_.size(); ++i)
			{
				/*out_file << i << std::endl;*/
				curr_face_ = iter->second.faces_.at(i);
				face_sqr_ = curr_face_->square();
				curr_normal_ = curr_face_->normal_vector();

				reconstructor_->reconstruct_vals(curr_face_, rvars_);
				for (size_t j = 0; j < 2; ++j)
				{
					vols_[j] = curr_face_->neighbour(j);
					vars_[j] = vols_[j]->variable();
					v_n_[j] = 0.0;
					v_sqr_ = 0.0;
					pressure_[j] = rvars_[j]->pressure();
					temperature_[j] = rvars_[j]->temperature();
					density_[j] = rvars_[j]->density();

					for (size_t k = 0; k < dim_; ++k)
					{
						vel_[j][k] = rvars_[j]->velocity(k);
						v_n_[j] += vel_[j][k] * curr_normal_->component(k);
						v_sqr_ += vel_[j][k] * vel_[j][k];
					}

					a_[j] = std::sqrt(gamma * rm * temperature_[j]);
					H_[j] = cp * temperature_[j] + 0.5 * v_sqr_;
					M_[j] = v_n_[j] / a_[j];
				}

				if (M_[1] >= 1.0)
				{
					M_sign_[1] = M_[1];
					p_sign_[1] = pressure_[1];
				}
				else if (M_[1] <= -1.0)
				{
					M_sign_[1] = 0.0;
					p_sign_[1] = 0.0;
				}
				else 
				{
					M_sign_[1] = 0.25 * (M_[1] + 1.0) * (M_[1] + 1.0);
					p_sign_[1] = 0.25 * pressure_[1] * (M_[1] + 1.0) * (M_[1] + 1.0) * (2.0 - M_[1]);
				}

				if (M_[0] >= 1.0)
				{
					M_sign_[0] = 0.0;
					p_sign_[0] = 0.0;
				}
				else if (M_[0] <= -1.0)
				{
					M_sign_[0] = M_[0];
					p_sign_[0] = pressure_[0];
				}
				else
				{
					M_sign_[0] = -0.25 * (M_[0] - 1.0) * (M_[0] - 1.0);
					p_sign_[0] = 0.25 * pressure_[0] * (M_[0] - 1.0) * (M_[0] - 1.0) * (2.0 + M_[0]);
				}

				M_f_ = M_sign_[0] + M_sign_[1];
				p_f_ = p_sign_[0] + p_sign_[1];
				
				if (M_f_ > 0.0)
					sign = 1;
				else
					sign = 0;

				temp = M_f_ * a_[sign] * density_[sign];
				for (size_t j = 0; j < dim_; ++j)
				{
					p_flow_[j] = p_f_ * curr_normal_->component(j);
					flow_[j] =  temp * vel_[sign][j] + p_flow_[j];
				}
				flow_[dim_] = temp;
				flow_[dim_plus] = temp * H_[sign];
				/*out_file << iter->first << " " << i << " " << H_Roe << " " << a_Roe << " " << v_Roe_normal << " "<< v_Roe_sqr << " ";*/
				sign = 1;
				for (size_t j = 0; j < 2; ++j)
				{
					sign = (-1) * sign;
					for (size_t k = 0; k < w_size; ++k)
					{
						vars_[j]->residual(k) += sign * flow_[k] * face_sqr_;
						/*out_file << sign << " " << flow[k] << " " << face_sqr << " " << vars[j]->residual(k) << std::endl;*/
					}
					//out_file << sign << " " << H[j] << a[j] << normal_velocity[j] << " " << v_sqr[j] << " ";
				}

				//for (size_t j = 0; j < w_size; ++j)
				//	out_file << flow[j] << " ";
				//out_file << std::endl;
			}
		}
	}
}

void psstConvectiveAUSM::get_data(psstGridFV & grid)
{
	grid_ = &grid;
}
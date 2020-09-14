#include "psstFiniteVolumeGrid.h"
#include "psstBoundaryCondition.h"
#include "psstConvectiveFlow.h"
#include "psstGradCalculation.h"
#include "psstLocalTimeCalculation.h"
#include "psstRungeKutta.h"
#include "psstScalarLimiters.h"
#include "psstValueReconstructor.h"
#include "psstVariable.h"
#include "psstOutputManager.h"
#include "psstFiniteVolumeGridTools.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <memory>
#include <vector>
#include <chrono>

double get_parameter(std::ifstream & in_file)
{
	double x = 0.0;
	in_file >> x;
	while (in_file.get() != '\n' && !in_file.eof())
		continue;
	return x;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstGridFV																																															|
|																																																		|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void psstGridFV::solver_parameters()
{
	std::string file_name;
	file_name = settings_file_name_ + ".inp";
	std::ifstream in_file(file_name);
	if (!in_file.is_open())
		throw std::exception("Error: Could not found setting file!");

	size_t init_mode = static_cast<size_t>(get_parameter(in_file));
	private_initialize_field(init_mode);
	private_load_boundary();

	double CFL = get_parameter(in_file);
	std::cout << "CFL number: " << CFL << std::endl;
	size_t lc = static_cast<size_t>(get_parameter(in_file));
	if (lc == 0)
	{
		if (std::abs(ref_vel_) < 1.0e-10)
			ref_vel_ = 1.0;
		time_step_calculator_ = std::make_shared<psstUniformTime>(ref_vel_);
		time_step_calculator_->CFL() = CFL;
		time_step_calculator_->get_data(*this);
		time_step_calculator_->initialize();
		std::cout << "Local Time Calculator: psstUniformTime" << std::endl;
	}
	else if (lc == 1)
	{
		time_step_calculator_ = std::make_shared<psstLocalAccTime1>();
		time_step_calculator_->CFL() = CFL;
		time_step_calculator_->get_data(*this);
		time_step_calculator_->initialize();
		std::cout << "Locat Time Calculator: psstLocalAccTime1" << std::endl;
	}
	else if (lc == 2)
	{
		time_step_calculator_ = std::make_shared<psstLocalAccTime2>();
		time_step_calculator_->CFL() = CFL;
		time_step_calculator_->get_data(*this);
		time_step_calculator_->initialize();
		std::cout << "Locat Time Calculator: psstLocalAccTime2" << std::endl;
	}

	std::cout << "Start from " << iter_ << " iteration" << std::endl;
	size_t max_iter = static_cast<size_t>(get_parameter(in_file));
	max_iter_ = max_iter;
	max_iter_ += iter_ - 1;
	std::cout << "Max number of iteration: " << max_iter_ << std::endl;

	size_t iter_step = static_cast<size_t>(get_parameter(in_file));
	iter_step_ = iter_step;
	std::cout << "Save iteration step " << iter_step_ << std::endl;

	size_t rk_level = static_cast<size_t>(get_parameter(in_file));
	std::cout << "Runge Kutta level: " << rk_level << std::endl;
	
	size_t rk_mode = static_cast<size_t>(get_parameter(in_file));
	if (rk_mode == 0)
	{
		std::cout << "Runge Kutta Standart Coefficients" << std::endl;
		RK_.standart_coefficients(rk_level);
		std::string temp_line;
		std::getline(in_file, temp_line);
	}
	else
	{
		std::cout << "Runge Kutta Custom Coefficients" << std::endl;
		std::vector<double> rk_coeff(rk_level);
		for (size_t i = 0; i < rk_level - 1; ++i)
		{
			in_file >> rk_coeff[i];
		}
		rk_coeff[rk_level - 1] = get_parameter(in_file);
		RK_.custom_coefficients(rk_coeff);
	}
	std::cout << "Runge Kutta coefficients:" << std::endl;
	for (size_t i = 0; i < RK_.level(); ++i)
	{
		std::cout << RK_[i] << " ";
	}

	size_t mode = static_cast<size_t>(get_parameter(in_file));
	size_t order = static_cast<size_t>(get_parameter(in_file));
	std::shared_ptr<psstInterfaceValueReconstructor> reconstructor;
	if (order == 1)
	{
		reconstructor = std::make_shared<psstFirstOrderReconstructor>();
		std::cout << "Scheme order: first order" << std::endl;
		calculate_grad_ = false;
	}
	else if (order == 2)
	{
		reconstructor = std::make_shared<psstReconstructorDM>(rm_, dimension_);
		std::cout << "Scheme order: second order" << std::endl;
		std::cout << "Using DM reconstruction" << std::endl;
		calculate_grad_ = true;
	}
	else
	{
		std::exception("Error! Unknown scheme order!");
	}
	std::cout << std::endl;

	switch (mode)
	{
		case 1:
			convective_flow_ = std::make_shared<psstConvectiveAUSM>();
			std::cout << "Using AUSM convective flow scheme" << std::endl;
		 break;
		case 2:
			convective_flow_ = std::make_shared<psstConvectiveHLL>();
			std::cout << "Using HLL convective flow scheme" << std::endl;
			break;
		case 3:
			convective_flow_ = std::make_shared<psstConvectiveHLLC>();
			std::cout << "Using HLLC convective flow scheme" << std::endl;
			break;
		default:
			throw std::exception("Error: Unknow convective flow scheme!");
	}
	convective_flow_->get_data(*this);

	size_t grad = static_cast<size_t>(get_parameter(in_file));
	if (grad == 1)
	{
		grad_calculator_ = std::make_shared <psstGreenGauss>();
		std::cout << "Using Green Gauss gradient method" << std::endl;
	}
	else
	{
		grad_calculator_ = std::make_shared<psstLeastSquare>();
		std::cout << "Using Least Square gradient method" << std::endl;
	}
	grad_calculator_->get_data(*this);
	grad_calculator_->initialize();

	size_t limiterType = static_cast<size_t>(get_parameter(in_file));
	if (limiterType == 1)
	{
		std::shared_ptr<psstInterfaceScalarLimiter> limiter = std::make_shared<psstLimiterVanAlbada>();
		reconstructor->attach_limiter(limiter);
		std::cout << "Using limiter VanAlbada" << std::endl;
	}
	else if (limiterType == 2)
	{
		std::shared_ptr<psstInterfaceScalarLimiter> limiter = std::make_shared<psstLimiterMinMod>();
		reconstructor->attach_limiter(limiter);
		std::cout << "Using limiter MinMod" << std::endl;
	}
	else if (limiterType == 3)
	{
		std::shared_ptr<psstInterfaceScalarLimiter> limiter = std::make_shared<psstLimiterVanLeer>();
		reconstructor->attach_limiter(limiter);
		std::cout << "Using limiter VanLeer" << std::endl;
	}
	
	convective_flow_->attach_reconstructor(reconstructor);
	convective_flow_->initialize();
	in_file.close();
}

std::vector<std::shared_ptr<psstInterfaceVertex>> * psstGridFV::vertex_storage()
{
	std::vector<std::shared_ptr<psstInterfaceVertex>> * vs_ptr = &vertex_storage_;
	return vs_ptr;
}

const std::vector<std::shared_ptr<psstInterfaceVertex>> * psstGridFV::vertex_storage() const
{
	const std::vector<std::shared_ptr<psstInterfaceVertex>> * vs_ptr = &vertex_storage_;
	return vs_ptr;
}

std::map<size_t, psstFaceStorage> * psstGridFV::face_storage()
{
	std::map<size_t, psstFaceStorage> * fs_ptr = &face_storage_;
	return fs_ptr;
}

const std::map<size_t, psstFaceStorage> * psstGridFV::face_storage() const
{
	const std::map<size_t, psstFaceStorage> * fs_ptr = &face_storage_;
	return fs_ptr;
}

std::map<size_t, psstVolumeStorage> * psstGridFV::volume_storage()
{
	std::map<size_t, psstVolumeStorage> * cs_ptr = &volume_storage_;
	return cs_ptr;
}

const std::map<size_t, psstVolumeStorage> * psstGridFV::volume_storage() const
{
	const std::map<size_t, psstVolumeStorage> * cs_ptr = &volume_storage_;
	return cs_ptr;
}

void psstGridFV::attach_output_manager(std::shared_ptr<psstInterfaceOutputManager> output_manager)
{
	if (output_manager.get() == 0)
		throw std::invalid_argument("Error: Expected output manager but null_ptr found!");
	psstInterfaceOutputManager * manager = output_manager->clone();
	output_manager_.reset(manager);
}

void psstGridFV::save_data() const
{
	if (output_manager_.get() == 0)
		throw std::exception("Error: Output manager is not attached!");
	output_manager_->output_data(*this);
	private_save_solution();
}

void psstGridFV::explicit_solve()
{
	double E = 0.0;
	double v_sqr = 0.0;
	double temp = 0.0;
	double density = 0.0;

	
	std::string file_name = settings_file_name_ + "_residual.dat";
	std::ofstream out_file;
	if (iter_ == 1)
	{
		out_file.open(file_name);
	}
	else
	{
		out_file.open(file_name, std::ios::app);
	}
	out_file.precision(3);
	out_file.scientific;

	size_t dim = dimension();
	std::map<size_t, psstVolumeStorage>::iterator iter = volume_storage_.begin();
	std::shared_ptr<psstInterfaceVolume> curr_vol = 0;
	std::shared_ptr<psstInterfaceVariable> curr_var = 0;
	std::shared_ptr<psstInterfaceBoundaryCondition> dyn_cond = 0;
	std::shared_ptr<psstInterfaceBoundaryCondition> t_cond = 0;
	size_t dim_plus = dim + 1;
	size_t w_size = dim + 2;

	auto begin = std::chrono::steady_clock::now();
	for (iter; iter != volume_storage_.end(); ++iter)
	{
		for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
		{
			curr_vol = iter->second.volumes_.at(i);
			curr_var = curr_vol->variable();
			v_sqr = 0.0;
			density = curr_var->density();
			for (size_t j = 0; j < dim; ++j)
			{
				temp = curr_var->velocity(j);
				curr_var->w(j) = density * temp;
				curr_var->residual(j) = 0.0;
				v_sqr += temp * temp;
			}
			E = cv_ * curr_var->temperature() + 0.5 * v_sqr;
			curr_var->residual(dim) = 0.0;
			curr_var->w(dim) = density;
			curr_var->residual(dim_plus) = 0.0;
			curr_var->w(dim_plus) = density * E;
		}
	}

	std::vector<double> res(w_size, 0.0);
	double vol = 0.0;
	double vel = 0.0;
	for (iter_; iter_ <= max_iter_; ++iter_)
	{
		if (calculate_grad_)
			grad_calculator_->calculate();
		for (size_t k = 0; k < RK_.level(); ++k)
		{
			convective_flow_->calculate();
			iter = volume_storage_.begin();
			for (iter; iter != volume_storage_.end(); ++iter)
			{
				if (iter->second.bc_type_id_ == 1)
				{
					for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
					{
						curr_vol = iter->second.volumes_.at(i);
						curr_var = curr_vol->variable();
						vol = curr_vol->volume();

						time_step_calculator_->calculate_at(curr_vol);
						for (size_t j = 0; j < w_size; ++j)
						{
							if (k == 0)
								curr_var->prev_w(j) = curr_var->w(j);
							curr_var->w(j) = curr_var->prev_w(j) -curr_var->time_step() * RK_[k] * curr_var->residual(j) / vol;

							/*curr_var->w(j) = curr_var->w(j) - curr_var->time_step() * RK_[k] * curr_var->residual(j) / vol;*/
							if (res[j] < std::abs(curr_var->residual(j)))
								res[j] = std::abs(curr_var->residual(j));
							curr_var->residual(j) = 0.0;
						}

						density = curr_var->w(dim);
						curr_var->density() = density;
						v_sqr = 0.0;
						for (size_t j = 0; j < dim; ++j)
						{
							curr_var->velocity(j) = curr_var->w(j) / density;
							vel = curr_var->velocity(j);
							v_sqr += vel * vel;
						}
						curr_var->temperature() = ((curr_var->w(dim_plus) / density) - (v_sqr*0.5)) / cv_;
						curr_var->pressure() = density * rm_ * curr_var->temperature();
					}
				}
			}

			iter = volume_storage_.begin();
			for (iter; iter != volume_storage_.end(); ++iter)
			{
				if (iter->second.bc_type_id_ == 1 || iter->second.bc_type_id_ == 32)
					continue;

				dyn_cond = iter->second.condition_->dynamic_cond;
				t_cond = iter->second.condition_->thermal_cond;
				for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
				{
					curr_var = iter->second.volumes_.at(i)->variable();
					dyn_cond->calculate_at(i);
					t_cond->calculate_at(i);
					curr_var->density() = curr_var->pressure() / (curr_var->temperature() * rm_);
				}
			}
		}

		std::cout << iter_ << " ";
		out_file << iter_ << " ";
		auto end = std::chrono::steady_clock::now();
		auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		auto time = elapsed_ms.count() / 1000.0;
		std::cout << std::setprecision(11) << std::scientific << time << " ";
		out_file << std::setprecision(11) << std::scientific << time << " ";

		for (size_t j = 0; j < w_size; ++j)
		{
			std::cout << std::setprecision(3) << std::scientific << " " << res[j];
			out_file << std::setprecision(3) << std::scientific << " "  << res[j];
			res[j] = 0.0;
		}
		std::cout << "\n";
		
		if ((iter_ % iter_step_) == 0)
		{
			std::cout << "Saving after " << iter_ << " iterations..." << std::endl;
			out_file << std::endl;
			save_data();
			//output_manager_->output_data(*this);
			//private_save_solution();
		}
		else
		{
			out_file << "\n";
		}
	}
	out_file.close();
	std::cout << "Final saving... " << std::endl;
	save_data();
}

void psstGridFV::private_load_boundary()
{
	std::string file_name;
	file_name = settings_file_name_ + ".bvp";
	std::ifstream in_file;
	in_file.open(file_name);
	if (!(in_file.is_open()))
		throw std::exception("Error! Could not open boundary condition file!");

	std::string zone_name;
	size_t bc_id;
	std::map<size_t, psstVolumeStorage>::iterator iter;
	std::cout << "Creating boundaries..." << std::endl;
	//Считывание динамических условий
	iter = volume_storage_.begin();
	while (!in_file.eof())
	{
		in_file >> zone_name >> bc_id;
		iter = volume_storage_.begin();
		while (iter->second.name_ != zone_name)
		{
			++iter;
			if (iter == volume_storage_.end())
				throw std::exception("Error: Could not found zone!");
		}
		if (iter->second.bc_type_id_ == 1 || iter->second.bc_type_id_ == 2 || iter->second.bc_type_id_ == 31 || iter->second.bc_type_id_ == 32)
			throw std::exception("Error: Could not attach boundary condition to inner zone!");

		std::shared_ptr<psstBoundary> b_cond = std::make_shared<psstBoundary>();
		enum dyn_bc_id
		{
			ss_inlet = 23,	//сверхзвуковое входное условие
			ss_outlet = 26,	//сверхзвуковое выходное условие
			wall_no_slip = 15,	//стенка с проскальзыванием
			inlet = 27,	//дозвуковое выходное условие
			outlet = 25	//дозвуковое выходное условие
		};

		if (bc_id == ss_inlet)
		{
			std::vector<double> velocity(dimension());
			double p;
			ref_vel_ = 0.0;
			for (size_t i = 0; i < dimension(); ++i)
			{
				in_file >> velocity[i];
				ref_vel_ += velocity[i] * velocity[i];
			}
			ref_vel_ = std::sqrt(ref_vel_);
			p = get_parameter(in_file);
			b_cond->dynamic_cond = std::make_shared<psstSupersonicInlet>(velocity,p);			
		}
		else if (bc_id == ss_outlet)
		{
			b_cond->dynamic_cond = std::make_shared<psstSupersonicOutlet>();
		}
		else if (bc_id == wall_no_slip)
		{
			b_cond->dynamic_cond = std::make_shared<psstInviscidWall>();
		}
		else if (bc_id == inlet)
		{
			std::vector<double> velocity(dimension());
			ref_vel_ = 0.0;
			for (size_t i = 0; i < dimension() - 1; ++i)
			{
				in_file >> velocity[i];
				ref_vel_ += velocity[i] * velocity[i];
			}
			velocity[dimension() - 1] = get_parameter(in_file);
			ref_vel_ += velocity[dimension() - 1] * velocity[dimension() - 1];
			ref_vel_ = std::sqrt(ref_vel_);
			b_cond->dynamic_cond = std::make_shared<psstInlet>(velocity);
		}
		else if (bc_id == outlet)
		{
			double p = 0.0;
			p = get_parameter(in_file);
			b_cond->dynamic_cond = std::make_shared<psstOutlet>(p);
		}
		b_cond->dynamic_cond->get_data(&volume_storage_[iter->first]);
		b_cond->dynamic_cond->initialize();
		iter->second.condition_ = b_cond;
	}
	in_file.close();

	//Считывание тепловых условий
	file_name = settings_file_name_ + ".btq";
	in_file.open(file_name);
	if (!(in_file.is_open()))
		throw std::exception("Error! Could not open boundary condition file!");
	
	while (!in_file.eof())
	{
		
		in_file >> zone_name >> bc_id;
		iter = volume_storage_.begin();
		while (iter->second.name_ != zone_name)
		{
			++iter;
			if (iter == volume_storage_.end())
				throw std::exception("Error: Could not found zone!");
		}
		if (iter->second.bc_type_id_ == 1 || iter->second.bc_type_id_ == 2 || iter->second.bc_type_id_ == 31 || iter->second.bc_type_id_ == 32)
			throw std::exception("Error: Could not attach boundary condition to inner zone!");

		enum term_bcID
		{
			wall_T = 11,	//стенка. T задается
			wall_Q = 18,	//стенка. Задается тепловой поток
			thermal_inlet = 27,	//входное условие
			thermal_outlet = 25	//выходное условие
		};
		std::shared_ptr<psstBoundary> b_cond;
		b_cond = iter->second.condition_;
		if (bc_id == wall_T)
		{
			double T = 0.0;
			T = get_parameter(in_file);
			double c = std::sqrt(gamma_ * rm_ * T);
			if (ref_vel_ < c)
				ref_vel_ = c;
			b_cond->thermal_cond = std::make_shared<psstThermalWall_T>(T);
		} else if (bc_id == wall_Q)
		{
			b_cond->thermal_cond = std::make_shared<psstThermalWall_Q>();
		}
		else if (bc_id == thermal_inlet)
		{
			double T = 0.0;
			T = get_parameter(in_file);
			double c = std::sqrt(gamma_ * rm_ * T);
			if (ref_vel_ < c)
				ref_vel_ = c;
			b_cond->thermal_cond = std::make_shared<psstThermalInlet>(T);
		}
		else if (bc_id == thermal_outlet)
		{
			b_cond->thermal_cond = std::make_shared<psstThermalOutlet>();
		}
		b_cond->thermal_cond->get_data(&volume_storage_[iter->first]);
		b_cond->thermal_cond->initialize();

		std::shared_ptr<psstInterfaceVariable> var = 0;
		for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
		{
			var = iter->second.volumes_.at(i)->variable();
			var->density() = var->pressure() / (rm_ * var->temperature());
		}
	}

	
	in_file.close();
	
}

void psstGridFV::private_load_parameter(std::ifstream & in_file)
{
	size_t zone_id = 0;
	size_t vol_size = 0;
	size_t dim_plus = dimension_ + 1;
	std::vector<std::shared_ptr<psstInterfaceVolume>> * vol = 0;
	std::shared_ptr<psstInterfaceVariable> var = 0;
	std::vector<double> v_init(dimension_, 0.0);
	double p_init = 0.0;
	double T_init = 0.0;
	double density = 0.0;

	in_file >> iter_;
	while (!in_file.eof())
	{
		in_file >> zone_id;
		if (zone_id == 0)
			break;

		vol = &volume_storage_.at(zone_id).volumes_;
		in_file >> vol_size;
		if (vol_size != vol->size())
			throw std::exception("Error! Different size of storage found!");

		for (size_t i = 0; i < vol_size; ++i)
		{
			var = vol->at(i)->variable();

			in_file >> p_init;
			var->pressure() = p_init;

			in_file >> T_init;
			var->temperature() = T_init;

			density = p_init / (rm_ * T_init);
			var->density() = density;

			for (size_t j = 0; j < dimension_; ++j)
			{
				in_file >> v_init[j];
				var->velocity(j) = v_init[j];
			}
		}

	}
}

void psstGridFV::private_save_solution() const
{
	std::ofstream out_file;
	std::string file_name = settings_file_name_ + ".sav";
	out_file.open(file_name);
	if (!out_file.is_open())
		throw std::exception("Error: Could not found *.sav file!");
	out_file.precision(11);

	out_file << iter_ << std::endl;
	std::map<size_t, psstVolumeStorage>::const_iterator iterator;
	iterator = volume_storage_.cbegin();
	size_t w_size = dimension_ + 2;
	for (iterator; iterator != volume_storage_.cend(); ++iterator)
	{
		if (iterator->second.bc_type_id_ != 1)
			continue;

		out_file << iterator->first << std::endl;
		out_file << iterator->second.volumes_.size() << std::endl;
		std::shared_ptr<psstInterfaceVariable> var = 0;
		for (size_t i = 0; i < iterator->second.volumes_.size(); ++i)
		{
			var = iterator->second.volumes_[i]->variable();
			out_file << std::fixed << var->pressure() << " ";
			out_file << std::fixed << var->temperature() << " ";
			for (size_t j = 0; j < dimension_; ++j)
				out_file << std::fixed << var->velocity(j) << " ";
			out_file << std::endl;
		}
	}
	out_file << 0;
	out_file.close();
}

void psstGridFV::private_initialize_field(size_t mode)
{
	std::string file_name;
	file_name = settings_file_name_ + ".prop";

	std::ifstream in_file(file_name);
	if (!in_file.is_open())
		throw std::exception("Error: Could not open fluid properties file");
	molecular_mass_ = get_parameter(in_file);
	cp_ = get_parameter(in_file);
	double R = 8.31;
	rm_ = R / molecular_mass_;
	gamma_ = cp_ / (cp_ - rm_);
	cv_ = cp_ / gamma_;
	in_file.close();


	if (mode == 0)
	{	
		size_t dim = dimension();
		size_t dim_plus = dim + 1;
		iter_ = 1;
		file_name = settings_file_name_ + ".par";
		in_file.open(file_name);
		if (!(in_file.is_open()))
			throw std::exception("Error! Could not open initial parameters file!");
		std::cout << "Init field..." << std::endl;
		std::vector<double> vel_init(dimension());
		double T_init = 0.0;
		double p_init = 0.0;
		double density = 0.0;
		for (size_t j = 0; j < dimension() - 1; ++j)
		{
			in_file >> vel_init[j];
		}
		vel_init[dimension() - 1] = get_parameter(in_file);
		T_init = get_parameter(in_file);
		p_init = get_parameter(in_file);
		in_file.close();
		
		std::shared_ptr<psstInterfaceVariable> var = 0;
		std::map<size_t, psstVolumeStorage>::iterator iter = volume_storage_.begin();
		for (iter; iter != volume_storage_.end(); ++iter)
		{
			if (iter->second.bc_type_id_ == 1)
			{
				for (size_t i = 0; i < iter->second.volumes_.size(); ++i)
				{
					var = iter->second.volumes_.at(i)->variable();
					
					for (size_t j = 0; j < dim; ++j)
					{
						var->velocity(j) = vel_init[j];
					}
					var->temperature() = T_init;
					var->pressure() = p_init;
					density = p_init / (rm_ * T_init);
					var->density() = density;
				}
			}
		}
	}
	else if (mode == 1)
	{
		file_name = settings_file_name_ + ".sav";
		in_file.open(file_name);
		if (!(in_file.is_open()))
			throw std::exception("Error! Could not open initial parameters file!");
		std::cout << "Loading field..." << std::endl;
		private_load_parameter(in_file);
		in_file.close();
	}
	else
	{
		throw std::exception("Error! Unknown command detected!");
	}
}
#ifndef _PSST_FINITE_VOLUME_GRID_H_
#define _PSST_FINITE_VOLUME_GRID_H_

#include <string>
#include <memory>
#include <vector>
#include <map>

#include "psstFiniteVolumeGridTools.h"
#include "psstRungeKutta.h"

class psstInterfaceGradientCalc;
class psstInterfaceLocalTime;
class psstInterfaceConvectiveFlow;
class psstInterfaceOutputManager;
class psstGridFV
{
private:
	size_t dimension_;
	std::vector<std::shared_ptr<psstInterfaceVertex>> vertex_storage_;
	std::map<size_t, psstFaceStorage> face_storage_;
	size_t face_number_;
	std::map<size_t, psstVolumeStorage> volume_storage_;
	size_t inner_volume_number_;
	
	double molecular_mass_;
	double min_vol_;
	double gamma_;
	double rm_;
	double cp_;
	double cv_;
	double ref_vel_;

	bool calculate_grad_;
	size_t iter_;
	size_t max_iter_;
	size_t iter_step_;

	psstRungeKutta RK_;
	std::shared_ptr<psstInterfaceGradientCalc> grad_calculator_;
	std::shared_ptr<psstInterfaceLocalTime> time_step_calculator_;
	std::shared_ptr<psstInterfaceConvectiveFlow> convective_flow_;
	std::shared_ptr<psstInterfaceOutputManager> output_manager_;
	
	std::string settings_file_name_;
	
	void private_load_boundary();
	void private_load_parameter(std::ifstream & inFile);
	void private_save_solution() const;
	void private_initialize_field(size_t mode);
public:
	virtual void settings_file_name(const std::string & file_name) { settings_file_name_ = file_name; };
	virtual void solver_parameters();
	virtual std::vector<std::shared_ptr<psstInterfaceVertex>> * vertex_storage();
	virtual const std::vector<std::shared_ptr<psstInterfaceVertex>> * vertex_storage() const;
	
	virtual size_t number_of_faces() const { return face_number_; };
	virtual void set_number_of_faces(size_t number_of_face) { face_number_ = number_of_face; };
	virtual std::map<size_t, psstFaceStorage> * face_storage();
	virtual const std::map<size_t, psstFaceStorage> * face_storage() const;
	virtual size_t number_of_inner_volumes() const { return inner_volume_number_; };
	virtual void set_number_of_inner_volumes(size_t number_of_volumes) { inner_volume_number_ = number_of_volumes; };

	virtual std::map<size_t, psstVolumeStorage> * volume_storage();
	virtual const std::map<size_t, psstVolumeStorage> * volume_storage() const;

	virtual void attach_output_manager(std::shared_ptr<psstInterfaceOutputManager> output_manager);
	virtual void save_data() const;

	virtual size_t dimension() const { return dimension_; };
	virtual void set_dimension(size_t dim) { dimension_ = dim; };
	
	virtual void set_min_vol(double vol) { min_vol_ = vol; };
	virtual double min_vol() const { return min_vol_;};
	
	virtual void set_gamma(double gamma) {gamma_ = gamma; };
	virtual double gamma() const { return gamma_; };

	virtual void set_molecular_mass(double molecular_mass) { molecular_mass_ = molecular_mass; };
	virtual double molecular_mass() const { return molecular_mass_; };

	virtual double rm() const { return rm_; };
	virtual double cp() const { return cp_; };
	virtual double cv() const { return cv_; };
	virtual void explicit_solve();
	virtual ~psstGridFV() {};
};

#endif

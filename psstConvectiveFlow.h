#ifndef _PSST_CONVECTIVE_FLOW_H_
#define _PSST_CONVECTIVE_FLOW_H_

#include <memory>
#include <vector>

double average_by_Roe(double left_val, double right_val, double left_density, double right_density);
double average_by_Roe(double left_density, double right_density);

class psstInterfaceValueReconstructor;
class psstGridFV;
class psstInterfaceVariable;
class psstInterfaceVolume;
class psstInterfaceVector;
class psstInterfaceFace;
class psstInterfaceConvectiveFlow
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceConvectiveFlow * clone() const = 0;
	virtual void attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor> ) = 0;
	virtual void calculate() = 0;
	virtual void get_data(psstGridFV & grid) = 0;
	virtual ~psstInterfaceConvectiveFlow() {};
};

class psstConvectiveHLL : public psstInterfaceConvectiveFlow
{
private:
	std::shared_ptr<psstInterfaceValueReconstructor> reconstructor_;
	psstGridFV * grid_;

	double face_sqr_;
	double pressure_[2];
	double temperature_[2];
	double density_[2];
	std::vector<std::vector<double>> vel_;
	std::vector<double> v_Roe_;
	double v_n_[2];
	double H_[2];
	double E_[2];
	double a_[2];
	double s_[2];
	double a_Roe_;
	double H_Roe_;
	double v_Roe_n_;
	double v_sqr_;
	size_t dim_;
	std::shared_ptr<psstInterfaceFace> curr_face_;
	std::shared_ptr<const psstInterfaceVector> curr_normal_;
	std::shared_ptr<psstInterfaceVolume> vols_[2];
	std::shared_ptr<psstInterfaceVariable> vars_[2];
	std::vector<std::shared_ptr<psstInterfaceVariable>> rvars_;

	std::vector<std::vector<double>> flow_i_;
	std::vector<std::vector<double>> w_i_;
	std::vector<double> flow_;
public:
	virtual void initialize();
	virtual psstInterfaceConvectiveFlow * clone() const;
	virtual void attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor>);
	virtual void calculate();
	virtual void get_data(psstGridFV & grid);
	virtual ~psstConvectiveHLL() {};
};

class psstConvectiveHLLC : public psstInterfaceConvectiveFlow
{
private:
	std::shared_ptr<psstInterfaceValueReconstructor> reconstructor_;
	psstGridFV * grid_;

	double face_sqr_;
	double pressure_[2];
	double temperature_[2];
	double density_[2];
	std::vector<std::vector<double>> vel_;
	std::vector<double> v_Roe_;
	double v_n_[2];
	double H_[2];
	double E_[2];
	double a_[2];
	double s_[2];
	double a_Roe_;
	double H_Roe_;
	double v_Roe_n_;
	double v_sqr_;
	double s_hllc_;
	size_t dim_;
	std::shared_ptr<psstInterfaceFace> curr_face_;
	std::shared_ptr<const psstInterfaceVector> curr_normal_;
	std::shared_ptr<psstInterfaceVolume> vols_[2];
	std::shared_ptr<psstInterfaceVariable> vars_[2];
	std::vector<std::shared_ptr<psstInterfaceVariable>> rvars_;

	std::vector<std::vector<double>> flow_i_;
	std::vector<std::vector<double>> w_i_;
	std::vector<double> flow_;
public:
	virtual void initialize();
	virtual psstInterfaceConvectiveFlow * clone() const;
	virtual void attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor>);
	virtual void calculate();
	virtual void get_data(psstGridFV & grid);
	virtual ~psstConvectiveHLLC() {};
};

class psstConvectiveAUSM : public psstInterfaceConvectiveFlow
{
private:
	std::shared_ptr<psstInterfaceValueReconstructor> reconstructor_;
	psstGridFV * grid_;

	double face_sqr_;
	double pressure_[2];
	double temperature_[2];
	double density_[2];
	std::vector<std::vector<double>> vel_;
	double v_n_[2];
	double M_[2];
	double H_[2];
	double a_[2];
	double M_sign_[2];
	double p_sign_[2];
	size_t dim_;
	double v_sqr_;
	double p_f_;
	double M_f_;

	std::shared_ptr<psstInterfaceFace> curr_face_;
	std::shared_ptr<const psstInterfaceVector> curr_normal_;
	std::shared_ptr<psstInterfaceVolume> vols_[2];
	std::shared_ptr<psstInterfaceVariable> vars_[2];
	std::vector<std::shared_ptr<psstInterfaceVariable>> rvars_;
	std::vector<double> p_flow_;
	std::vector<double> flow_;
public:
	virtual void initialize();
	virtual psstInterfaceConvectiveFlow * clone() const;
	virtual void attach_reconstructor(std::shared_ptr<psstInterfaceValueReconstructor>);
	virtual void calculate();
	virtual void get_data(psstGridFV & grid);
	virtual ~psstConvectiveAUSM() {};
};

#endif
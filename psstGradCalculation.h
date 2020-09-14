#ifndef _PSST_GRAD_CALCULATION_H_
#define _PSST_GRAD_CALCULATION_H_
#include <string>
#include <memory>
#include <vector>
#include <map>

class psstGridFV;
class psstInterfaceFace;
class psstInterfaceVolume;
class psstInterfaceVector;
class psstInterfaceVariable;
class psstInterfaceVertex;
class psstInterfaceGradientCalc
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceGradientCalc * clone() const = 0;
	virtual void get_data(psstGridFV & grid) = 0;
	virtual void calculate() const = 0;
	virtual ~psstInterfaceGradientCalc() {};
};

class psstGreenGauss : public psstInterfaceGradientCalc
{
private:
	mutable std::shared_ptr<psstInterfaceFace> curr_face_;
	mutable std::shared_ptr<const psstInterfaceVector> normal_;
	mutable std::shared_ptr<psstInterfaceVariable> var_[2];
	mutable std::shared_ptr<psstInterfaceVolume> vol_[2];
	size_t dim_;
	mutable double geom_par_[2];
	mutable std::vector<std::vector<double *>> ptr_vars_;
	mutable	std::vector<std::vector<std::shared_ptr<psstInterfaceVector>>> ptr_grads_;
	mutable std::vector<double> temp_vals_;
	std::vector<std::vector<double>> distances_;
	std::vector<double> beta_;
	psstGridFV * grid_ = 0;
public:
	virtual void initialize();
	virtual psstInterfaceGradientCalc * clone() const;
	virtual void get_data(psstGridFV & grid);
	virtual void calculate() const ;
	virtual ~psstGreenGauss() {};
};

class psstLeastSquare : public psstInterfaceGradientCalc
{
private:
	psstGridFV * grid_ = 0;

	std::vector<std::vector<double>> r_;
	std::vector<double> beta_;
	mutable std::vector < std::shared_ptr<psstInterfaceVector>> ptr_grads_;
	enum index
	{
		a = 0,
		b = 1,
		d = 2,
		c = 3,
		e = 4,
		g = 5,
	};
public:
	virtual void initialize();
	virtual psstInterfaceGradientCalc * clone() const;
	virtual void get_data(psstGridFV & grid);
	virtual void calculate() const;
	virtual ~psstLeastSquare() {};
};

#endif
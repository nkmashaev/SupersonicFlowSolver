#include "psstFiniteVolumeGridTools.h"
#ifndef _PSST_VARIABLE_H_
#define _PSST_VARIABLE_H_

class psstInterfaceVector;
class psstInterfaceVariable
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceVariable * clone() const = 0;
	virtual void copy(const psstInterfaceVariable &) = 0;
	virtual size_t dimension() const = 0;

	virtual double & velocity(size_t i) = 0;
	virtual const double & velocity(size_t i) const = 0;
	virtual double & temperature() = 0;
	virtual const double & temperature() const = 0;
	virtual double & pressure() = 0;
	virtual const double & pressure() const = 0;
	virtual double & density() = 0;
	virtual const double & density() const = 0;
	
	virtual std::shared_ptr<psstInterfaceVector> & grad_velocity(size_t i) = 0;
	virtual std::shared_ptr<const psstInterfaceVector> grad_velocity(size_t i) const = 0;
	virtual std::shared_ptr<psstInterfaceVector> & grad_temperature() = 0;
	virtual std::shared_ptr<const psstInterfaceVector> grad_temperature() const = 0;
	virtual std::shared_ptr<psstInterfaceVector> & grad_pressure() = 0;
	virtual std::shared_ptr<const psstInterfaceVector> grad_pressure() const = 0;
	virtual std::shared_ptr<psstInterfaceVector> & grad_density() = 0;
	virtual std::shared_ptr<const psstInterfaceVector> grad_density() const = 0;
	virtual void nullify_grad() = 0;

	virtual double & w(size_t i) = 0;
	virtual const double & w(size_t i) const = 0;
	
	virtual double & residual(size_t i) = 0;
	virtual const double & residual(size_t i) const = 0;

	virtual double & time_step() = 0;
	virtual const double & time_step() const = 0;

	virtual double & prev_w(size_t i) = 0;
	virtual const double & prev_w(size_t i) const = 0;
	virtual ~psstInterfaceVariable() {};
};

class psstVariable : public psstInterfaceVariable
{
private:
	std::vector<double> velocity_;
	double temperature_;
	double pressure_;
	double density_;
	std::vector<std::shared_ptr<psstInterfaceVector>> grad_velocity_;
	std::shared_ptr<psstInterfaceVector> grad_temperature_;
	std::shared_ptr<psstInterfaceVector> grad_pressure_;
	std::shared_ptr<psstInterfaceVector> grad_density_;
	std::vector<double> w_;
	std::vector<double> prev_w_;
	std::vector<double> residual_;
	double time_step_;
	size_t dimension_;

	void private_set_dimension(size_t);
public:
	psstVariable();
	psstVariable(size_t dim);
	psstVariable(const psstVariable &);
	virtual void initialize();
	virtual psstInterfaceVariable * clone() const;
	virtual void copy(const psstInterfaceVariable &);
	virtual size_t dimension() const { return dimension_; };

	virtual double & velocity(size_t i);
	virtual const double & velocity(size_t i) const;
	virtual double & temperature();
	virtual const double & temperature() const;
	virtual double & pressure();
	virtual const double & pressure() const;
	virtual double & density();
	virtual const double & density() const;

	virtual std::shared_ptr<psstInterfaceVector> & grad_velocity(size_t i);
	virtual std::shared_ptr<const psstInterfaceVector> grad_velocity(size_t i) const;
	virtual std::shared_ptr<psstInterfaceVector> & grad_temperature();
	virtual std::shared_ptr<const psstInterfaceVector> grad_temperature() const;
	virtual std::shared_ptr<psstInterfaceVector> & grad_pressure();
	virtual std::shared_ptr<const psstInterfaceVector> grad_pressure() const;
	virtual std::shared_ptr<psstInterfaceVector> & grad_density();
	virtual std::shared_ptr<const psstInterfaceVector> grad_density() const;
	virtual void nullify_grad();

	virtual double & w(size_t i);
	virtual const double & w(size_t i) const;

	virtual double & residual(size_t i);
	virtual const double & residual(size_t i) const;

	virtual double & time_step();
	virtual const double & time_step() const;

	virtual double & prev_w(size_t i) { return prev_w_.at(i); };
	virtual const double & prev_w(size_t i) const { return prev_w_.at(i); };

	psstVariable & operator=(const psstVariable & var) { copy(var); return *this; };
	virtual ~psstVariable() ;
};


#endif

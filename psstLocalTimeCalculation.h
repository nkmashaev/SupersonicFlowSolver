#ifndef _PSST_LOCAL_TIME_CALCULATOR_H_
#define _PSST_LOCAL_TIME_CALCULATOR_H_

#include <vector>
#include <memory>

class psstGridFV;
class psstInterfaceVolume;
class psstInterfaceFace;
class psstInterfaceVariable;
class psstInterfaceVector;
class psstInterfaceLocalTime
{
public:
	virtual void initialize() = 0;
	virtual psstInterfaceLocalTime * clone() const = 0;

	virtual double & CFL() = 0;
	virtual const double & CFL() const = 0;
	virtual void get_data(psstGridFV & grid) = 0;
	virtual void calculate_at(std::shared_ptr<psstInterfaceVolume>) const = 0;
	virtual ~psstInterfaceLocalTime() {};
};
	
class psstUniformTime : public psstInterfaceLocalTime
{
private:
	double ref_velocity_;
	double CFL_;
	psstGridFV * grid_ = 0;

	psstUniformTime();
public:
	psstUniformTime(double velocity);
	virtual void initialize();
	virtual psstInterfaceLocalTime * clone() const;

	virtual double & CFL() { return CFL_; };
	virtual const double & CFL() const { return CFL_; };
	virtual void get_data(psstGridFV & grid);
	virtual void calculate_at(std::shared_ptr<psstInterfaceVolume>) const {};
	virtual ~psstUniformTime() {};
};

struct psstVolumeStorage;
class psstLocalAccTime1 : public psstInterfaceLocalTime
{
private:
	double CFL_;
	psstGridFV * grid_ = 0;

	double gamma_;
	double rm_;
	size_t dim_;

	mutable std::shared_ptr<const psstInterfaceVariable> neighbour_var_ = 0;
	mutable std::shared_ptr<psstInterfaceVariable> var_ = 0;
	mutable std::shared_ptr<psstInterfaceFace> curr_face_ = 0;
	mutable std::shared_ptr<const psstInterfaceVector> normal_ = 0;
public:
	virtual void initialize();
	virtual psstInterfaceLocalTime * clone() const;

	virtual double & CFL() { return CFL_; };
	virtual const double & CFL() const { return CFL_; };
	virtual void get_data(psstGridFV & grid);
	virtual void calculate_at(std::shared_ptr<psstInterfaceVolume>) const;
	virtual ~psstLocalAccTime1() {};
};

class psstLocalAccTime2 : public psstInterfaceLocalTime
{
private:
	double CFL_;
	psstGridFV * grid_ = 0;

	double gamma_;
	double rm_;
	size_t dim_;

	mutable std::vector<double> s_;
	mutable std::shared_ptr<psstInterfaceVariable> var_ = 0;
	mutable std::shared_ptr<psstInterfaceFace> curr_face_ = 0;
	mutable std::shared_ptr<const psstInterfaceVector> normal_ = 0;
public:
	virtual void initialize();
	virtual psstInterfaceLocalTime * clone() const;

	virtual double & CFL() { return CFL_; };
	virtual const double & CFL() const { return CFL_; };
	virtual void get_data(psstGridFV & grid);
	virtual void calculate_at(std::shared_ptr<psstInterfaceVolume>) const;
	virtual ~psstLocalAccTime2() {};
};
#endif

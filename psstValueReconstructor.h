#ifndef _PSST_VOLUME_RECONSTRUCTOR_H_
#define _PSST_VOLUME_RECONSTRUCTOR_H_

#include "psstFiniteVolumeGridTools.h";

#include <memory>
#include <vector>

class psstInterfaceScalarLimiter;
class psstInterfaceVariable;
class psstInterfaceValueReconstructor
{
public:
	virtual psstInterfaceValueReconstructor * clone() const = 0;
	virtual void attach_limiter(std::shared_ptr<const psstInterfaceScalarLimiter>) = 0;
	virtual void reconstruct_vals(std::shared_ptr<const psstInterfaceFace>, std::vector<std::shared_ptr<psstInterfaceVariable>> & rvars) = 0;
	virtual ~psstInterfaceValueReconstructor() {};
};

class psstFirstOrderReconstructor : public psstInterfaceValueReconstructor
{
private:
	std::shared_ptr<const psstInterfaceScalarLimiter> limiter_;
public:
	virtual psstInterfaceValueReconstructor * clone() const;
	virtual void attach_limiter(std::shared_ptr<const psstInterfaceScalarLimiter>);
	virtual void reconstruct_vals(std::shared_ptr<const psstInterfaceFace>, std::vector<std::shared_ptr<psstInterfaceVariable>> & rvars);
	virtual ~psstFirstOrderReconstructor() {};
};

class psstReconstructorDM : public psstInterfaceValueReconstructor
{
private:
	std::shared_ptr<const psstInterfaceScalarLimiter> limiter_;
	double rm_;
	size_t dim_;
	std::vector<double> R_LR_;
	std::vector<std::vector<double *>> ptr_vars_;
	std::vector<std::vector<std::shared_ptr<const psstInterfaceVector>>> ptr_grads_;
	psstReconstructorDM();
public:
	psstReconstructorDM(double rm, size_t dim) : 
		rm_(rm), 
		dim_(dim),
		limiter_(0)
		{ 
			R_LR_.resize(dim);
			ptr_vars_ = std::vector<std::vector<double *>>(2, std::vector<double*>(dim + 3, 0));
			ptr_grads_ = std::vector<std::vector<std::shared_ptr<const psstInterfaceVector>>>(2, std::vector<std::shared_ptr<const psstInterfaceVector>>(dim + 3, 0));
		};
	virtual psstInterfaceValueReconstructor * clone() const;
	virtual void attach_limiter(std::shared_ptr<const psstInterfaceScalarLimiter>);
	virtual void reconstruct_vals(std::shared_ptr<const psstInterfaceFace>, std::vector<std::shared_ptr<psstInterfaceVariable>> & rvars);
	virtual ~psstReconstructorDM() {};
};

#endif
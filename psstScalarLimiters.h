#ifndef _PSST_SCALAR_LIMITERS_H_
#define _PSST_SCALAR_LIMITERS_H_

class psstInterfaceScalarLimiter
{
public:
	virtual psstInterfaceScalarLimiter * clone() const = 0;
	virtual double get_value(double numerator, double denumerator) const = 0;
	virtual ~psstInterfaceScalarLimiter() {};
};

class psstLimiterVanAlbada : public psstInterfaceScalarLimiter
{
public:
	virtual psstInterfaceScalarLimiter  * clone() const;
	virtual double get_value(double numerator, double denumerator) const;
	virtual ~psstLimiterVanAlbada() {};
};

class psstLimiterMinMod : public psstInterfaceScalarLimiter
{
public:
	virtual psstInterfaceScalarLimiter  * clone() const;
	virtual double get_value(double numerator, double denumerator) const;
	virtual ~psstLimiterMinMod() {};
};

class psstLimiterVanLeer : public psstInterfaceScalarLimiter
{
public:
	virtual psstInterfaceScalarLimiter  * clone() const;
	virtual double get_value(double numerator, double denumerator) const;
	virtual ~psstLimiterVanLeer() {};
};
#endif

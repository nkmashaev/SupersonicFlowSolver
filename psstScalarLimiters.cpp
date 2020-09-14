#include "psstScalarLimiters.h"
#include <cmath>

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceScalarLimiter																																											|
|	psstLimiterVanAlbada																																												|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstInterfaceScalarLimiter  * psstLimiterVanAlbada::clone() const
{
	psstInterfaceScalarLimiter  * new_limiter = new psstLimiterVanAlbada(*this);
	return new_limiter;
}

double psstLimiterVanAlbada::get_value(double numerator, double denumerator) const
{
	double eps = 1.0e-10;
	if ((numerator > 0.0 && denumerator < 0.0) || (numerator < 0.0 && denumerator > 0.0))
		return 0.0;

	if (std::abs(denumerator) < eps)
		return 1.0;

	double r = 0.0;
	double res = 0.0;
	r = numerator / denumerator;
	res = (r * r + r) / (r * r + 1.0);
	return res;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceScalarLimiter																																											|
|	psstLimiterMinMod																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstInterfaceScalarLimiter  * psstLimiterMinMod::clone() const
{
	psstInterfaceScalarLimiter  * new_limiter = new psstLimiterMinMod(*this);
	return new_limiter;
}

double psstLimiterMinMod::get_value(double numerator, double denumerator) const
{
	double eps = 1.0e-10;
	if ((numerator > 0.0 && denumerator < 0.0) || (numerator < 0.0 && denumerator > 0.0))
		return 0.0;

	if (std::abs(denumerator) < eps)
		return 1.0;

	double r = 0.0;
	double res = 0.0;
	r = numerator / denumerator;
	res = r;
	if (r > 1.0)
		res = 1.0;
	return res;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|	psstInterfaceScalarLimiter																																											|
|	psstLimiterVanLeer																																													|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

psstInterfaceScalarLimiter  * psstLimiterVanLeer::clone() const
{
	psstInterfaceScalarLimiter  * new_limiter = new psstLimiterVanLeer(*this);
	return new_limiter;
}


double psstLimiterVanLeer::get_value(double numerator, double denumerator) const
{
	double eps = 1.0e-10;
	if ((numerator > 0.0 && denumerator < 0.0) || (numerator < 0.0 && denumerator > 0.0))
		return 0.0;

	if (std::abs(denumerator) < eps)
		return 2.0;

	double r = 0.0;
	double res = 0.0;
	r = numerator / denumerator;

	res = 2.0 * r / (1.0 + r);
	return res;
}

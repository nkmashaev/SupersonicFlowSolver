#ifndef _PSST_RungeKutta_H_
#define _PSST_RungeKutta_H_

#include<vector>

class psstRungeKutta
{
private:
	size_t rk_lvl_;
	std::vector<double> rk_coeff_;
public:
	psstRungeKutta() :
		rk_lvl_(1),
		rk_coeff_(1)
	{
		standart_coefficients(level());
	}

	psstRungeKutta(size_t lvl) :
		rk_lvl_(lvl),
		rk_coeff_()
	{
		standart_coefficients(level());
	}


	void standart_coefficients(size_t lvl);
	void custom_coefficients(const std::vector<double> & coeff);
	size_t level() const { return rk_lvl_; };

	const double & operator[](size_t index);
	virtual ~psstRungeKutta() {};
};

#endif

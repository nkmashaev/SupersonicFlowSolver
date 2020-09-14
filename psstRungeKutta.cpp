#include "psstRungeKutta.h"

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|																																																		|
|	psstRungeKutta																																														|																																																		|
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

void  psstRungeKutta::standart_coefficients(size_t level)
{
	rk_coeff_.clear();
	switch (level)
	{
	case 1:
		rk_coeff_.resize(level);
		rk_coeff_[0] = 1.0;
		break;
	case 2:
		rk_coeff_.resize(level);
		rk_coeff_[0] = 0.4242;
		rk_coeff_[1] = 1.0;
		break;
	case 3:
		rk_coeff_.resize(level);
		rk_coeff_[0] = 0.1918;
		rk_coeff_[1] = 0.4929;
		rk_coeff_[2] = 1.0000;
		break;
	case 4:
		rk_coeff_.resize(level);
		rk_coeff_[0] = 0.1084;
		rk_coeff_[1] = 0.2602;
		rk_coeff_[2] = 0.5052;
		rk_coeff_[3] = 1.0000;
		break;
	case 5:
		rk_coeff_.resize(level);
		rk_coeff_[0] = 0.0695;
		rk_coeff_[1] = 0.1602;
		rk_coeff_[2] = 0.2898;
		rk_coeff_[3] = 0.5060;
		rk_coeff_[4] = 1.0000;
		break;
	default:
		throw std::exception("Error: Unknow RK order!");
	}
	rk_lvl_ = level;
}

void psstRungeKutta::custom_coefficients(const std::vector<double> & coeff)
{
	if (coeff.size() == 0)
		throw std::exception("Error! Expected Runge Kutta coefficients but empty vector found!");
	rk_lvl_ = coeff.size();
	rk_coeff_ = coeff;
}

const double & psstRungeKutta::operator[](size_t index)
{
	return rk_coeff_.at(index);
}
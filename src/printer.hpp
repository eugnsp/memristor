#pragma once
#include <iomanip>
#include <ostream>

class Printer
{
public:
	Printer(std::ostream& os) : os_(os)
	{
		os_ << std::setprecision(4);
	}

	void num(unsigned int n)
	{
		if (n % 20 == 0)
			header();

		os_ << std::left << std::setw(4) << n + 1 << ' ' << std::flush;
	}

	template<typename T>
	void value(T value, T scale = 1, unsigned int width = 10)
	{
		os_ << std::setw(width - 1) << value / scale << ' ' << std::flush;
	}

	void endl()
	{
		os_ << std::endl;
	}

private:
	void header()
	{
		os_ << "-----------------------------------------------------------------------------------------------------------\n"
			<< "#    Time      Fil.vol   Ext.bias  Syst.bias Current   Syst.res. Lim.res.  MC est.      MC act.   Num.vacs.\n"
			<< "     [msec]    [arb.]    [V]       [V]       [μA]      [kΩ]      [kΩ]      [μsec/step]  [steps]            \n"
			<< "-----------------------------------------------------------------------------------------------------------"
			<< std::endl;
	}

private:
	std::ostream& os_;
};


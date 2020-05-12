/*
 * hydrate.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_HYDRATE_HH_
#define PARAMETERS_HYDRATE_HH_

class Hydrate
{
private:

public:

	constexpr static double grainDensity = 920.;

double density( ) const {
		/* unit -> kg/m^3 */

		double rho_h = grainDensity ;

		return rho_h;

	}

	double molarMass() const
	{
		/* unit -> kg/mol */
		return 119.5/1000;
	}

	double hydrationNumber() const
	{
		return 5.90; //6.176; used for kiel triax test
	}

	double thermalConductivity( double T, double Ppore ) const
	{
		double kth;
		/* mu: unit -> W.m^-1 K^-1 */

		kth = 0.5 ;

		return kth;
	}

	double Cp( double T, double Ppore ) const
	{
		double Cp;
		/* mu: unit -> J/kg.K */

		Cp = 2216.0;

		return Cp;
	}

	double Cv( double T, double Ppore ) const
	{
		double Cv;
		/* mu: unit -> J/kg.K */

		Cv = Cp( T, Ppore ) ;

		return Cv;
	}

};



#endif /* PARAMETERS_HYDRATE_HH_ */

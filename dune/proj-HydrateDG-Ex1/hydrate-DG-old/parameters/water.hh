/*
 * water.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_WATER_HH_
#define PARAMETERS_WATER_HH_

class Water
{

public:

	double density( double T, double Pw ) const {
		double rho;
		/* rho: unit -> kg/m^3 */

		rho=1000.0;

		return rho;
	}

	double vaporDensity( double T, double Pg ) const {
		double vapRho;
		/* rho: unit -> kg/m^3 */

		vapRho=0.0022*Pg/T;
		// REFERENCE : http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html
		// Here Pg is taken instead of water vapor partial pressure ( i.e. Y_H2O * Pg ), because, in here it assumes Y_H2O = 1
		// Then, to get actual vapRho while calculating avg Gas Density, we multiply this vapRho with Y_H2O
		// This is basically done so that density is called only as a function of T and Pg

		return vapRho;
	}

	double dynamicViscosity( double T, double Pw ) const {
		double mu;
		/* mu: unit -> Pa.s */

//		mu = 0.5 * 1.0e-3 ;

		// REFERENCE:
		double mu_0 = 0.001792 ; // kg/m/s
		double a = - 1.94 ;
		double b = - 4.80 ;
		double c =  6.74 ;
		double T0 = 273.15 ; // K
		double Tr = T0/T ;
		mu = mu_0 * exp( a + b * Tr + c * Tr*Tr );

		return mu;
	}

	double molarMass( ) const {
		/* unit -> kg/mol */
		return 18.0/1000;
	}

	double thermalConductivity( double T, double Pw ) const {
		// ref:http://people.ucsc.edu/~bkdaniel/WaterProperties.html
		// making a logarithmic fit to this data
		/* 0 C   -> 0.561
		 * 20 C  -> 0.5984
		 * 40 C  -> 0.6305
		 * 60 C  -> 0.6543
		 * 80 C  -> 0.670
		 * 100 C -> 0.6791
		 */

		double kth;
		/* mu: unit -> W.m^-1 K^-1 */

		kth = 0.3834 * log(T) - 1.581 ;

		return kth;
	}

	double Cp( double T, double Pw ) const {
		double Cp;
		/* mu: unit -> J*kg^-1*K^-1 */

		Cp = 4186.0 ;

		return Cp;
	}

	double Cv( double T, double Pw ) const {
		double Cv;
		/* mu: unit -> J*kg^-1*K^-1 */

		Cv = Cp( T, Pw );

		return Cv;
	}

	double saturatedVaporPressure( double T ) const {
		double A = 0.;
		double B = 0.;
		double C = 0.;

		T = T - 273.15;
		if ( T <= 100.0 ){
			A = 8.07131;
			B = 1730.63;
			C = 233.426;
		}
		else if ( T > 100.0 ){
			A = 8.14019;
			B = 1810.94;
			C = 244.485;
		}
		double Psat = pow( 10, ( A - B/( C + T ) ) ); /* [mmHg] */

		/*  760 mmHg = 101.325 kPa */
		return Psat * (101.325*1000.0/760.0) ;	/* [Pa] */
	}

};


#endif /* PARAMETERS_WATER_HH_ */

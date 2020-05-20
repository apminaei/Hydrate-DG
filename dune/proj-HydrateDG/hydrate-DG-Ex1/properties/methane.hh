/*
 * methane.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_METHANE_HH_
#define PARAMETERS_METHANE_HH_


/*___________________________________________________________
 *
 * CHECK OUT PROPERTIES ON :-
 * http://www.ddbst.de/en/EED/PCP/PCPindex.php#Methane
 *___________________________________________________________
 */

class Methane
{
private:
	/* http://en.wikipedia.org/wiki/Gas_constant */
	constexpr static double Ru = 8.3144; /* [J*mol^-1*K^-1] */
	Indices indices;
public:

	double molarMass() const {
		return 16.04 * 1.0e-3; 	/* [kg/mol] */
	}

	double accentricityFactor()const {
		return 0.011 	 ;
	}

	double criticalTemperature( ) const {
		return -82.7 + 273.15 ; /* [K] */
	}

	double criticalPressure( ) const {
		return 45.96 * 1.0e5 ; /* [Pa] */
	}

	double density(double T, double Pg, double z_CH4) const {
		// double rho;
		// double R_CH4 = Ru/molarMass();
		// /* rho: unit -> kg/m^3 */

		// rho = Pg/( z_CH4 * R_CH4 * T);

		return 19.605;
	}

	double densityAtSTP( ) const {
		double rho;
		double R_CH4 = Ru/molarMass();
		double std_T = 273.15 + 0.;
		double std_P = 1.013 * 1.e5;
		double z_CH4 = 1.0;

		rho = std_P/( z_CH4 * R_CH4 * std_T);

		return rho;
	}

	double dynamicViscosity(double T, double Pg) const {
		// double mu;
		// /* mu: unit -> Pa.s */
		// // Sutherland Correlation:
		// // ref: http://portal.tpu.ru/SHARED/n/NATASHA/Material/Tab3/Glava_1.pdf
		// double C = 162; // empirical constant
		// double mu_0 = 10.4 * 1.e-6 ; // Pa.s
		// // ref for mu_0 :http://www.pipeflowcalculations.com/tables/gas.php
		// mu = mu_0 * (273.15 + C ) * ( pow( (T/273.15), 1.5) / ( T + C ) ) ;

		return 1.1045e-5;
	}

	double thermalConductivity( double T, double Pg) const {
		// double kth;

		// // REFERENCE: " Thermal Conductivity of Methane for temperatures between 110 K and 310 K with Pressures upto 70 MPa
		// //			  author : H.M. Roder
		// //			  Journal: Journal of Thermophysics, volume 6, No 2, pages 119-142
		// // assumption: dilute gas , therefore, density effects are neglected
		// double A0 = - 0.8863333440 * 1.e-2 ;
		// double A1 =   0.2419639784 * 1.e-3 ;
		// double A2 = - 0.6997019196 * 1.e-6 ;
		// double A3 =   0.1224609018 * 1.e-8 ;

		// kth = A0 + A1 * T + A2 * T*T + A3 * T*T*T ; /* [W*m^-1*K^-1] */

		return 0.03107;
	}

	double Cp_ideal( double T, double Pg ) const {
		/* REF: 1D Modelling of Hydrate Decomposition in Porous Media, by F. Esmailzadeh, M.E. Zeighami, J. Fathi */
		// double Cp_i;

		// double A = 1.238 ;
		// double B = 0.00313 ;
		// double C = 7.905*1.0e-7 ;
		// double D = -6.858*1.0e-10 ;

		// Cp_i = A + B*T + C*T*T +D*T*T*T ;		/* [kJ/(kg*K)] */
		// Cp_i*1000.0
		return  2165.24;		/* [J/(kg*K)] */
	}

	double Cp_res( double T, double Pg, double z_CH4 ) const {

		// Based on Peng Robinson's EoS
		// REFERENCE:
		double omega = accentricityFactor();
		double Tc	 = criticalTemperature();
		double Pc 	 = criticalPressure();

		double kappa = 0.;
		if( omega <= 0.49 ){
			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
		}
		else{
			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
		}

		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );

		double b = 0.07780 * ( Ru * Tc / Pc );
		double a_T = 0.45724 * ( ( Ru * Ru * Tc * Tc ) / Pc ) * ac ;

		double da_T = kappa * ac * ( ( kappa / Tc ) - (( 1 + kappa )/( sqrt( T*Tc ) ) ) );
		double dda_T = ( kappa * ac * ( 1. + kappa ) ) / ( 2. * sqrt( T*Tc ));

		double A = ( a_T * Pg ) / ( pow( ( Ru * T ) , 2 ) ) ;
		double B = ( b * Pg ) / ( Ru * T ) ;
		double M = ( z_CH4*z_CH4 + 2*B*z_CH4 - B*B ) / ( z_CH4 - B ) ;
		double N = ( da_T * B ) / ( b * Ru );

		double Cp_res =   dda_T * ( T / (2*sqrt(2) * b ) ) * log((z_CH4+(sqrt(2)+1)*B)/(z_CH4-(sqrt(2)-1)*B))
						+ ( Ru * pow( M-N , 2 ) ) / ( M*M - 2.*A * (z_CH4+B) )
						- Ru ;

		return Cp_res;
	}

	double Cp( double T, double Pg, double z_CH4 ) const {
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double Cp;
		Cp = Cp_ideal( T, Pg ) + Cp_res( T, Pg, z_CH4 );
		return Cp;
	}

	double Cv(double T, double Pg, double z_CH4 ) const {
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double omega = accentricityFactor();
		double Tc	 = criticalTemperature();
		double Pc 	 = criticalPressure();

		double kappa = 0.;
		if( omega <= 0.49 ){
			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
		}
		else{
			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
		}

		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );

		double b = 0.07780 * ( Ru * Tc / Pc );
		double a_T = 0.45724 * ( ( Ru * Ru * Tc * Tc ) / Pc ) * ac ;

		double da_T = kappa * ac * ( ( kappa / Tc ) - (( 1 + kappa )/( sqrt( T*Tc ) ) ) );
		double dda_T = ( kappa * ac * ( 1. + kappa ) ) / ( 2. * sqrt( T*Tc ));

		double A = ( a_T * Pg ) / ( pow( ( Ru * T ) , 2 ) ) ;
		double B = ( b * Pg ) / ( Ru * T ) ;
		double M = ( z_CH4*z_CH4 + 2*B*z_CH4 - B*B ) / ( z_CH4 - B ) ;
		double N = ( da_T * B ) / ( b * Ru );

		double Cv =   Cp( T, Pg, z_CH4 ) - ( Ru * pow( M-N , 2 ) ) / ( M*M - 2.*A * (z_CH4+B) ) ;

		return Cv; /* [J/(kg*K)] */
	}

	double henrysConstant( double T ) const {
		/* Reference : Lide and Frederkise, 1995 */

		/* Henry's Law constant for solubility in water at 298.15 K */
		// double k0_H = 0.0014;	/* [mol*kg^-1*Bar^-1] */
		// /* Temperature dependant constant */
		// double dlnkH_d1byT = 1600.0;	/* [K] */

		// /* ref: webbook.nist.gov/cgi/cbook.cgi?1D=C74828&Mask=10 */
		// double kHenry = k0_H*exp( dlnkH_d1byT *( 1.0/T - 1.0/298.15 ) ); /* [mol*kg^-1*Bar^-1] */

		// return kHenry * molarMass() * 1.0e-5; /* [Pa^-1] */
		return 3.e-10;
	}

};

/*
 -  Critical point
 --------------------
 *	Critical temperature  : -82.7 °C
 *	Critical pressure  : 45.96 bar
 *
 -	Gaseous phase
 ---------------------
 *	Gas density (1.013 bar at boiling point) : 1.819 kg/m3
 *	Gas density (1.013 bar and 15 °C (59 °F)) : 0.68 kg/m3
 *	Compressibility Factor (Z) (1.013 bar and 15 °C (59 °F)) : 0.998
 *	Specific gravity (air = 1) (1.013 bar and 21 °C (70 °F)) : 0.55
 *	Specific volume (1.013 bar and 21 °C (70 °F)) : 1.48 m3/kg
 *	Heat capacity at constant pressure (Cp) (1 bar and 25 °C (77 °F)) : 0.035 kJ/(mol.K)
 *	Heat capacity at constant volume (Cv) (1 bar and 25 °C (77 °F)) : 0.027 kJ/(mol.K)
 *	Ratio of specific heats (Gamma:Cp/Cv) (1 bar and 25 °C (77 °F)) : 1.305454
 *	Viscosity (1.013 bar and 0 °C (32 °F)) : 0.0001027 Poise
 *	Thermal conductivity (1.013 bar and 0 °C (32 °F)) : 32.81 mW/(m.K)

*/


#endif /* PARAMETERS_METHANE_HH_ */

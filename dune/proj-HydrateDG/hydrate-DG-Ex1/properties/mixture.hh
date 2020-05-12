/*
 * mixture.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_MIXTURE_HH_
#define PARAMETERS_MIXTURE_HH_


class Mixture
{
private:
	Methane methane;
	Water water;

public:

	/* CALCULATION OF COMPOSITIONS */

	double fnc_CH4( double T, double Pg ) const {

		return 1.0/(Pg * methane.henrysConstant( T ));

	}

	double fnc_H2O( double T, double Pg ) const {

		return water.saturatedVaporPressure( T )/Pg;

	}

	/* MOLE FRACTIONS ( x -> liquid ; y -> gas ) */

	double mole_x_CH4( double T, double Pg ) const {

		return ( 1.0 - fnc_H2O( T,Pg ) )/( fnc_CH4( T,Pg ) - fnc_H2O( T,Pg ) );

	}

	double mole_x_H2O( double T, double Pg ) const {

		return 1.0 - mole_x_CH4( T, Pg ) ;//( 1.0 - fnc_CH4( T,Pg ) )/( fnc_H2O( T,Pg ) - fnc_CH4( T,Pg ) );

	}

	double mole_y_CH4( double T, double Pg ) const {

		return fnc_CH4( T,Pg )*( 1.0 - fnc_H2O( T,Pg ) )/( fnc_CH4( T,Pg ) - fnc_H2O( T,Pg ) );

	}

	double mole_y_H2O( double T, double Pg ) const {

		return 1.0 - mole_y_CH4( T, Pg ) ;//fnc_H2O( T,Pg )*( 1.0 - fnc_CH4( T,Pg ) )/( fnc_H2O( T,Pg ) - fnc_CH4( T,Pg ) );

	}


	/* MASS FRACTIONS ( X -> liquid ; Y -> gas ) */

	double X_CH4( double T, double Pg ) const {
		// M_i*x_i/(sum(M_i*xi))
		return methane.molarMass( ) * mole_x_CH4( T,Pg )/( methane.molarMass( ) * mole_x_CH4( T,Pg ) + water.molarMass( ) * mole_x_H2O( T,Pg ) );

	}

	double X_H2O( double T, double Pg ) const {
		// M_i*x_i/(sum(M_i*xi))
		return water.molarMass( ) * mole_x_H2O( T,Pg )/( methane.molarMass( ) * mole_x_CH4( T,Pg ) + water.molarMass( ) * mole_x_H2O( T,Pg ) );

	}

	double Y_CH4( double T, double Pg ) const {
		// M_i*x_i/(sum(M_i*xi))
		return methane.molarMass( ) * mole_y_CH4( T,Pg )/( methane.molarMass( ) * mole_y_CH4( T,Pg ) + water.molarMass( ) * mole_y_H2O( T,Pg ) );

	}

	double Y_H2O( double T, double Pg ) const {
		// M_i*x_i/(sum(M_i*xi))
		return methane.molarMass( ) * mole_y_H2O( T,Pg )/( methane.molarMass( ) * mole_y_CH4( T,Pg ) + water.molarMass( ) * mole_y_H2O( T,Pg ) );

	}

	/* MASS TRANSFER COEFFICIENTS */

	double binaryDiffCoeffInGas( double T, double Pg ) const {

		// double a0 = 0. ;
		// double a1 = 2.26e-9;
		// double a2 = 0.002554;

		// return ( a0 + a1*T + a2/Pg ) ;	/* m^2/s */
		return 0.637e-6;
	}

	double binaryDiffCoeffInLiquid( double T, double Pg ) const {

		// double A = 0.003475 ; 	/* K */
		// double B = 1.57e-5;		/* cm^2/s */

		// return B * exp(-A/T) * 1.0e-6 ;	/* m^2/s */
		return 1.57e-11;
	}

	double DiffCoeffSaltInLiquid( double T/*K*/, double Pg/*Pa*/ ) const {

		double D;
		/* m^2/s */

		double a0 = 0. ;
		double a1 = 2.26e-9;
		double a2 = 0.002554;

		D = 1.0e-9 ;

		return D; /* m^2/s */

	}


	/* AVERAGE MIXTURE PROPERTIES FOR EACH PHASE */

	double avgGasDensity( double T, double Pg, double Pw, double z_CH4 ) const {

		double rhoW = water.vaporDensity( T , Pg );
		double rhoNW = methane.density( T , Pg , z_CH4 );

		double avgRhoNW = rhoNW * Y_CH4( T,Pg ) + rhoW * Y_H2O( T,Pg ) ;

		return avgRhoNW;

	}

	double avgLiquidDensity( double T, double Pg, double Pw, double z_CH4 ) const {

		double rhoW = water.density( T,Pw );
		double rhoNW = methane.density( T,Pw,z_CH4 );

		double avgRhoW = rhoNW * X_CH4( T,Pg ) + rhoW * X_H2O( T,Pg ) ;

		return avgRhoW;

	}

	double avgGasCp( double T, double Pg, double Pw, double z_CH4 ) const {

		double CpW = water.Cp( T,Pg );
		double CpNW = methane.Cp( T,Pg,z_CH4 );

		double avgCpNW = Y_CH4( T,Pg ) * CpNW + Y_H2O( T,Pg ) * CpW ;

		return avgCpNW;

	}

	double avgLiquidCp( double T, double Pg, double Pw, double z_CH4 ) const {

		double CpW = water.Cp( T,Pw );
		double CpNW = methane.Cp( T,Pw,z_CH4 );

		double avgCpW = X_CH4( T,Pg ) * CpNW + X_H2O( T,Pg ) * CpW ;

		return avgCpW;

	}

	double avgGasCv( double T, double Pg, double Pw, double z_CH4 ) const {

		double CvW = water.Cv( T,Pg );
		double CvNW = methane.Cv( T,Pg,z_CH4 );

		double avgCvNW = Y_CH4( T,Pg ) * CvNW + Y_H2O( T,Pg ) * CvW ;

		return avgCvNW;

	}

	double avgLiquidCv( double T, double Pg, double Pw, double z_CH4 ) const {

		double CvW = water.Cp( T,Pw );
		double CvNW = methane.Cp( T,Pw,z_CH4 );

		double avgCvW = X_CH4( T,Pg ) * CvNW + X_H2O( T,Pg ) * CvW ;

		return avgCvW;

	}

	double avgGasViscosity( double T, double Pg, double Pw ) const {

		double muW = water.dynamicViscosity( T,Pg );
		double muNW = methane.dynamicViscosity( T,Pg );

		double avgMuNW = muNW * X_CH4( T,Pg ) + muW * X_H2O( T,Pg ) ;

		return avgMuNW;

	}

	double avgLiquidViscosity( double T, double Pg, double Pw ) const {

		double muW = water.dynamicViscosity( T,Pw );
		double muNW = methane.dynamicViscosity( T,Pw );

		double avgMuW = muNW * X_CH4( T,Pg ) + muW * X_H2O( T,Pg ) ;

		return avgMuW;

	}

};



#endif /* PARAMETERS_MIXTURE_HH_ */

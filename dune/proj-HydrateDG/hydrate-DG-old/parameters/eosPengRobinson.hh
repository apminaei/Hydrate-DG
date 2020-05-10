/*
 * eosPengRobinson.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_EOSPENGROBINSON_HH_
#define PARAMETERS_EOSPENGROBINSON_HH_


class PengRobinson{
private:
	constexpr static double R_u = 8.314462175;
	constexpr static double eps = 1.e-6;
	constexpr static int iterMax = 10;
	Methane methane;
public:

	void printName(){
		std::cout<<"EOS: PENG-ROBINSON " << std::endl;
	}

	std::vector<double> evaluateEoSParams( double T, double P )const{

		std::vector<double> PREoSParams(2);
		for(int i =0; i<PREoSParams.size();i++){
			PREoSParams[i] = 0.;
		}

		double omega = methane.accentricityFactor();
		double Tc	 = methane.criticalTemperature();
		double Pc 	 = methane.criticalPressure();

		double kappa = 0.;
		if( omega <= 0.49 ){
			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
		}
		else{
			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
		}

		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );

		double a_T = 0.45724 * ( ( R_u * R_u * Tc * Tc ) / Pc ) * ac ;
		double b = 0.07780 * ( R_u * Tc / Pc );

		PREoSParams[0] = a_T;
		PREoSParams[1] = b  ;

		return PREoSParams;

	}

	std::vector<double> polynomialCoefficients( double T, double P )const{

		std::vector<double> EoSParameters(2);
		for(int i =0; i<EoSParameters.size();i++){
			EoSParameters[i] = 0.;
		}
		EoSParameters = evaluateEoSParams( T,P );

		double a_T = EoSParameters[0];
		double b   = EoSParameters[1];

		double A = ( a_T * P ) / ( pow( ( R_u * T ) , 2 ) ) ;
		double B = ( b * P ) / ( R_u * T ) ;

		std::vector<double> Coeffs(4);
		for(int i =0; i<Coeffs.size();i++){
			Coeffs[i] = 0.;
		}
		// cubic equation: a0*z^3 + a1*z^2 + a2*z + a3

		Coeffs[0] = 1;
		Coeffs[1] = B - 1. ;
		Coeffs[2] = A - 2.*B - 3.*B*B ;
		Coeffs[3] = B*B*B + B*B - A*B ;

		return Coeffs;

	}

	double method( std::vector<double> Coeffs )const{

//		std::cout <<"**NEWTON RAPHSON METHOD CALLED FOR Z_CH4 CALCULATION **"<< std::endl;

		int polynomialOrder = Coeffs.size() - 1;
//		std::cout<< "PolynomialOrder = " << polynomialOrder << std::endl;

		double z = 1.;
		double z_up = 1.;
		double defect = 0.;

		int iter = 0;
		do{
			iter += 1;
			double f_z = 0.;
			double df_z = 0.;
			for( int i = 0 ; i<=polynomialOrder ; i++ ){
				f_z += Coeffs[i] * pow( z , polynomialOrder - i );
				df_z += ( polynomialOrder - i ) * Coeffs[i] * pow( z , (polynomialOrder - i) -1 ) ;
			}
			z_up = z - (f_z / df_z );
			defect = z_up - z;

//			std::cout<< "iteration number = " << iter << "  , defect = " << defect << "  , z = " << z_up << std::endl;

			z = z_up;
		}
		while( std::abs(defect)>eps and iter<iterMax );

//		std::cout <<"z_CH4 CALCULATED = " << z_up << std::endl;

		return z_up;
	}

	double evaluateCompressibilityFactor( double T, double P )const{

		std::vector<double> Coeffs(4);
		for(int i =0; i<Coeffs.size();i++){
			Coeffs[i] = 0.;
		}
		Coeffs = polynomialCoefficients(T,P);
		double compressibilityFactor = method( Coeffs );
		return compressibilityFactor ;
	}

};


#endif /* PARAMETERS_EOSPENGROBINSON_HH_ */

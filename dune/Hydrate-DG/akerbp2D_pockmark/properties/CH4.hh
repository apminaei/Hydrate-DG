/* ALL PARAMETERS ARE NONDIMENSIONAL */
class Methane
{
private:
	CharacteristicValues characteristicValue;

public:

	/* http://en.wikipedia.org/wiki/Gas_constant */
	constexpr static double Ru = 8.314462175; /* [J*mol^-1*K^-1] */

	double MolarMass() const {
		return 16.04 * 1.0e-3; 	/* [kg/mol] */
	}

	double AccentricityFactor()const {
		return 0.011 	 ;
	}

	double CriticalTemperature( ) const {
		return -82.7 + 273.15 ; /* [K] */
	}

	double CriticalPressure( ) const {
		return 45.96 * 1.0e5 ; /* [Pa] */
	}

	double Density(double T/*K*/, double Pg/*Pa*/, double z_CH4) const {

		double rho;
		/* rho: unit -> kg/m^3 */

		double R_CH4 = Ru/MolarMass();
		rho = Pg/( z_CH4 * R_CH4 * T);

		return 19.605/characteristicValue.density_c; /*ndim*/
	}

	double DynamicViscosity(double T/*K*/, double Pg/*Pa*/) const {

		double mu;
		/* mu: unit -> Pa.s */

		// Sutherland Correlation:
		// ref: http://portal.tpu.ru/SHARED/n/NATASHA/Material/Tab3/Glava_1.pdf
		double C = 162; // empirical constant
		double mu_0 = 1.0707e-5;
		// ref for mu_0 :http://www.pipeflowcalculations.com/tables/gas.php
		mu_0 *= (  1. - (1./(1.0707e-5)) * (  4.8134e-14 * Pg
											- 4.1719e-20 * Pg * Pg
											+ 7.3232e-28 * Pg * Pg * Pg )
				); // Pa.s -> ref: Frauenhofer Comsol Model
		mu = mu_0 * (273.15 + C ) * ( pow( (T/273.15), 1.5) / ( T + C ) ) ;

		return 1.1045e-5/characteristicValue.viscosity_c;/*ndim*/

	}

	double ThermalConductivity( double T/*K*/, double Pg/*Pa*/) const {

		double kth; /* [W*m^-1*K^-1] */

#ifdef FRAUENHOFER_MODEL
		// REFERENCE: Frauenhofer Comsol Model
		kth = 43.82e-3 ; /* [W*m^-1*K^-1] */
#else
		// REFERENCE: " Thermal Conductivity of Methane for temperatures between 110 K and 310 K with Pressures upto 70 MPa
		//			  author : H.M. Roder
		//			  Journal: Journal of Thermophysics, volume 6, No 2, pages 119-142
		// assumption: dilute gas , therefore, density effects are neglected
		double A0 = - 0.8863333440 * 1.e-2 ;
		double A1 =   0.2419639784 * 1.e-3 ;
		double A2 = - 0.6997019196 * 1.e-6 ;
		double A3 =   0.1224609018 * 1.e-8 ;
		kth = A0 + A1 * T + A2 * T*T + A3 * T*T*T ;
#endif
		
		return 0.03107/characteristicValue.thermalconductivity_c;/*ndim*/
	}

	double Cp_ideal( double T/*K*/, double Pg/*Pa*/ ) const {

		double Cp_i;
		/* [J/(kg*K)] */

#ifdef FRAUENHOFER_MODEL
		// REFERENCE: Frauenhofer COMSOL model
		Cp_i = 3274.0; /* [J/(kg*K)] */
#else
		/* REF: 1D Modelling of Hydrate Decomposition in Porous Media, by F. Esmailzadeh, M.E. Zeighami, J. Fathi */
		double A = 1.238 ;
		double B = 0.00313 ;
		double C = 7.905*1.0e-7 ;
		double D = -6.858*1.0e-10 ;
		Cp_i = ( A + B*T + C*T*T +D*T*T*T ) * 1000.0 ; /* [J/(kg*K)] */
#endif

		return Cp_i ; /* [J/(kg*K)] */
	}

	double Cp_res( double T/*K*/, double Pg/*Pa*/, double z_CH4 ) const {

		double Cp_res;
		/* [J/(kg*K)] */

#ifdef FRAUENHOFER_MODEL
		Cp_res = 0.;
#else
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double omega = AccentricityFactor();
		double Tc	 = CriticalTemperature();
		double Pc 	 = CriticalPressure();

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

		Cp_res =   dda_T * ( T / (2*sqrt(2) * b ) ) * log((z_CH4+(sqrt(2)+1)*B)/(z_CH4-(sqrt(2)-1)*B))
				+ ( Ru * pow( M-N , 2 ) ) / ( M*M - 2.*A * (z_CH4+B) )
				- Ru ;
#endif

		return Cp_res; /*[J/(kg.K)]*/
	}

	double Cp( double T/*K*/, double Pg/*Pa*/, double z_CH4 ) const {
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double Cp;
		/* [J/(kg*K)] */

		Cp = Cp_ideal( T, Pg ) + Cp_res( T, Pg, z_CH4 );		/* [J/(kg*K)] */
		
		return 2165.24/characteristicValue.specificheat_c; /*ndim*/
	}

	double Cv(double T/*K*/, double Pg/*Pa*/, double z_CH4 ) const {

		double Cv;
		/* [J/(kg*K)] */

		Cv/*[J/(kg.K)]*/ =   Cp( T, Pg, z_CH4 )*characteristicValue.specificheat_c; // NOTE: Cases are checked in Cp_ideal function!

#ifdef FRAUENHOFER_MODEL
		Cv += (-1.) * Ru/MolarMass();		/* [J/(kg*K)] */
#else
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double omega = AccentricityFactor();
		double Tc	 = CriticalTemperature();
		double Pc 	 = CriticalPressure();

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

		double nonidealfactor = ( pow( M-N , 2 ) ) / ( M*M - 2.*A * (z_CH4+B) );

		Cv += (-1.) * Ru * nonidealfactor ;
#endif

		return Cv/characteristicValue.specificheat_c;/*ndim*/
	}

	double SolubilityCoefficient( double T/*K*/, double S ) const {

		double kHenry;
		/* [Pa^-1] */

		// REF: SUGAR TOOLBOX
		Water water;
		double Tc /*[K]*/   = water.CriticalTemperature();	//critical temperature of water
		double Pc /*[MPa]*/ = water.CriticalPressure()/1e6;	//critical pressure of water
		double tau = 1 - T/Tc;    							//dimensionless temperature
		double Ts = T/Tc;         							//reduced temperature
		//std::cout<< "Ts = " << Ts << " Pc = " << Pc << " S = " << S << std::endl;
		// vapor pressure of water
		double a1 = -7.85951783;
		double a2 =  1.84408259;
		double a3 = -11.7866497;
		double a4 =  22.6807411;
		double a5 = -15.9618719;
		double a6 =  1.80122502;
		double pvap /*[MPa]*/ = Pc * exp( Tc/T * ( a1*tau + a2*pow(tau,1.5) + a3*pow(tau,3) + a4*pow(tau,3.5) 
															+ a5*pow(tau,4) + a6*pow(tau,7.5) ) );
		//double pvap /*[MPa]*/ = water.SaturatedVaporPressure( T/*K*/, S );
		//std::cout<< "pvap = " << pvap << std::endl;
		// Henry constant
		double A = -11.0094;
		double B = 4.8362;
		double C = 12.5220;
		kHenry /*MPa*/ = exp( log(pvap) + A/Ts + B*pow((1 - Ts),0.355)/Ts + C*exp(1 - Ts)*pow(Ts,(-0.41)) );
		kHenry *= 1.e6; /* [Pa] */
		// std::cout << kHenry << std::endl;
		// exit(0);
		return 3.1104e9/characteristicValue.P_c;//kHenry/characteristicValue.P_c; /*ndim*/
	}

};

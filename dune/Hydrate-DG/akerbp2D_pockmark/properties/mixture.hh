/* ALL PARAMETERS ARE NONDIMENSIONAL */
class Mixture
{
private:
	Methane methane;
	Water water;
	Salt salt;
	CharacteristicValues X_c;

public:

	/* MOLE FRACTIONS ( X -> liquid ; Y -> gas ) */

	std::vector<double> EquilibriumMoleFractions( double T/*K*/, double Pg/*Pa*/, double Sg, double Sw , double Xc, double z )const{

		double Y_H2O, Y_CH4, X_H2O, X_CH4;

		// NOTE: it is not necessary to check state dependence for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/methane.MolarMass());
		double f_CH4 /*ndim*/ = z*(Pg/X_c.P_c)/methane.SolubilityCoefficient(T,S);
		double f_H2O /*ndim*/ = water.SaturatedVaporPressure( T,S )/(Pg/X_c.P_c);

		// Y_H2O = ((1.-Xc)-f_CH4)/(f_H2O-f_CH4);
		// Y_CH4 = 1.-Y_H2O;
		// X_H2O = Y_H2O * f_H2O;
		// X_CH4 = 1. - Xc - X_H2O;
		double zP = 1./(1. - z * water.SaturatedVaporPressure( T,S )/ methane.SolubilityCoefficient(T,S) );
		Y_CH4 = zP  * (1. - (1.-Xc) * water.SaturatedVaporPressure( T,S )/(Pg/X_c.P_c));
		X_H2O = zP * (1. - Xc  -  f_CH4);
		
		Y_H2O = f_H2O * X_H2O;
		X_CH4 = Y_CH4 * f_CH4;
		
			
		
		// X_CH4 = f_CH4*Y_CH4;// (z*water.SaturatedVaporPressure( T,S ))*(f_H2O+ Xc -1.)/(methane.SolubilityCoefficient(T,S)-(z*water.SaturatedVaporPressure( T,S )));
		// X_H2O = f_H2O*Y_H2O;//1.-X_CH4-Xc;

		std::vector<double> X(Indices::numOfComps,0.);
		X[Indices::compId_XCH4] = X_CH4;
		X[Indices::compId_XH2O] = X_H2O;
		X[Indices::compId_YCH4] = Y_CH4;
		X[Indices::compId_YH2O] = Y_H2O;

		return X;
	}

	double YCH4( double X_CH4, double T, double Pg, double Xc, double z )const{

		double Y_CH4;

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/methane.MolarMass());
		Y_CH4 = X_CH4 * methane.SolubilityCoefficient(T,S) / ( z * Pg/X_c.P_c ) ;

		return Y_CH4;
	}

	double XH2O( double Y_H2O, double T, double Pg, double Xc )const{

		double X_H2O;

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/methane.MolarMass());
		X_H2O = Y_H2O * (Pg/X_c.P_c) / water.SaturatedVaporPressure( T,S );

		return X_H2O;
	}

	/* MASS TRANSFER COEFFICIENTS */

	double DiffCoeffH2OInGas( double T/*K*/, double Pg/*Pa*/ ) const {

		double D;
		/* m^2/s */

		double a0 = 0. ;
		double a1 = 2.26e-9;
		double a2 = 0.002554;

		D = ( a0 + a1*T + a2/Pg ) ;

		return 0.637e-6/X_c.dispersivity_c; /*ndim*/

	}

	double DiffCoeffCH4InLiquid( double T/*K*/, double Pw/*Pa*/ ) const {

		double D; /* m^2/s */

		double A = 0.003475 ; 	/* K */
		double B = 1.57e-5;		/* cm^2/s */

		D = B * exp(-A/T) * 1.0e-6;

		return 1.57e-11/X_c.dispersivity_c; /*ndim*/
	}

};

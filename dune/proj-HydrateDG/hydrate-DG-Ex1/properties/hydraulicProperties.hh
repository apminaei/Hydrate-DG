// /*
//  * hydraulicProperties.hh
//  *
//  *  Created on: Sep 22, 2016
//  *      Author: shubhangi
//  */

#ifndef PARAMETERS_HYDRAULICPROPERTIES_HH_
#define PARAMETERS_HYDRAULICPROPERTIES_HH_
template<typename GV,typename PTree>
class HydraulicProperties
{
private:
	const GV& gv;
	const PTree& ptree;
	const static int dim = GV::dimension;
	constexpr static double m 		 	= 3.;
	constexpr static double beta	 	= 1.;
	const static int numOfParams  = 5;
	  const static int id_Pentry 	= 0;
	  const static int id_lambda 	= 1;
	  const static int id_Swr 		= 2;
	  const static int id_Sgr 		= 3;
	  const static int id_beta		= 4;
	Parameters<PTree> parameter;
	MeshParameters<PTree> mesh;
	//Soil<GV,PTree> soil;
	CharacteristicValues characteristicValue;
public:

	HydraulicProperties (const GV& gv_,const PTree& ptree_)
	: gv( gv_ ),
	ptree(ptree_),
	parameter(ptree_),
	mesh(ptree_)
	{}
/* BROOKS-COREY */

/*****************************************************************************/
	std::vector<double> BrooksCoreyParameters( ) const {
		
		auto prop_L = parameter.layer_properties();
		double Pentry = 50000. ;
		double lambda = 1.2 ;
		double Sgr = 0. ;
		double Swr = 0. ;
		double beta = 0. ;

		std::vector<double> BCParams (numOfParams,0.);
		BCParams[id_Pentry] = prop_L[0][2] ; /*Pa*/
		BCParams[id_lambda] = prop_L[0][3] ;
		BCParams[id_Swr]    = prop_L[0][4] ;
		BCParams[id_Sgr]    = prop_L[0][5] ;
		BCParams[id_beta]   = prop_L[0][6] ;

		return BCParams;
	}

/*****************************************************************************/
	double Pentry() const{
		double Pentry;
		Pentry = BrooksCoreyParameters()[0];
		return Pentry ;
	}

	double lambda() const {
		double lambda;
		lambda = BrooksCoreyParameters()[1];
		return lambda ;
	}

	double Swr() const{
		double Swr;
		Swr = BrooksCoreyParameters()[2];
		return Swr ;
	}

	double Sgr() const{
		double Sgr;
		Sgr = BrooksCoreyParameters()[3];
		return Sgr ;
	}
	double betaf() const{
		double beta;
		beta = BrooksCoreyParameters()[4];
		return beta ;
	}
/*****************************************************************************/

	double effectiveSw( double Sw, double Sh ) const
	{
		double Sw_max = 1. - Sh - Sgr();
		double Swe = ( Sw - Swr() )/( Sw_max - Swr() );

		return  Swe ;
	}
	

	double dSwe_dSw( double Sw, double Sh ) const {
		double Sw_max = 1. - Sh - Sgr();
		double dSwe =  1./( Sw_max - Swr() );

		return dSwe;
	}

	double dSwe_dSh( double Sw, double Sh ) const {
		double Sw_max = 1. - Sh - Sgr();
		double Swe = ( Sw - Swr() )/( Sw_max - Swr() );
		double dSwe = Swe /( Sw_max - Swr() );

		return dSwe;
	}

	/* SUCTION/CAPILLARY PRESSURE */

	double suctionPressure( double Sw, double Sh ) const {
		double eta = (1/lambda());
		double Swe = effectiveSw( Sw,Sh );
		double Pc = 0.;
		double a = 0.05 ;
		double b = 0.9 ;

		if( Swe > a /*and Swe < b*/ ){
			Pc = Pentry() * pow( Swe, -eta );
		}
		else if ( Swe <= a ){
			double Pc_a = Pentry() * pow( a, -eta );
			double dPc_a = dPc_dSwe( a,0. ) ;
			Pc = Pc_a + dPc_a * ( Swe - a );
		}
		else if ( Swe >= b and Swe < 1.0 - 1.e-2 ){
			double Pc_b = Pentry() * pow( b, -eta );
			double Pc_1 = 0.+ 1.e-3 ;
			double slope = ( Pc_b - Pc_1 ) / ( b - 1. ) ;
			double deltaSwe = Swe - 1.0 ;
			if( deltaSwe < 1.e-3 ){
				deltaSwe = 0.;
			}
			Pc = Pc_1 + slope * deltaSwe ;
		}
		else if( Swe >= 1.-1.e-2 ){
			Pc = 0. + 1.e-3;
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::suctionPressure( Sw,Sh ) "
					 << "  , Swe = " << Swe
					 << "  , Sw  = " << Sw
					 << "  , Sh  = " << Sh
					 << "  , Pc  = " << Pc << std::endl;
			exit(0);
		}

		if( Pc < -1e-3 ){
			std::cout<< " Pc is -ve " << std::endl;
			std::cout<< " ERROR in HydraulicProperties::suctionPressure( Sw,Sh )"
					 << "  , Swe = " << Swe
					 << "  , Sw  = " << Sw
					 << "  , Sh  = " << Sh
					 << "  , Pc  = " << Pc << std::endl;
			exit(0);
		}

		if( Pc < 0. + 1.e-3 ){
			Pc = 0.+ 1.e-3;
		}
		return Pc;
	}


	double dPc_dSwe( double Sw, double Sh ) const {

		double eta = (1/lambda());
		double Swe = effectiveSw( Sw,Sh );
		double dPc = 0.;
		double a = 0.05 ;

		if( Swe > a ){
			dPc = Pentry() * (-1./lambda()) * std::pow( effectiveSw(Sw,Sh) , -(1./lambda()) - 1. );
		}
		else if ( Swe <= a ){
			double dPc_a  = Pentry() * (-1./lambda()) * std::pow( a , -(1./lambda()) - 1. ) ;
			double ddPc_a = Pentry() * (-1./lambda()) * (-1./lambda()-1.) * std::pow( a , -(1./lambda()) - 2. );
			dPc = dPc_a + ddPc_a * ( Swe - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::dPc_dSwe( Sw,Sh ) "
					 << "  , Swe = " << Swe
					 << "  , Sw  = " << Sw
					 << "  , Sh  = " << Sh
					 << "  , dPc  = " << dPc << std::endl;
			exit(0);
		}

		return dPc;
	}

	double PcSF1( double Sh ) const
	{
		double PcSF=0. ;
		double eta = (lambda() * m - 1.)/(lambda()*m);
		double a = 0.95;
		if( Sh < a ){
			PcSF =  pow( 1.0 - Sh , -eta );
		}
		else if( Sh >= a ){
			double PcSF_a = std::pow( 1.0 - a , -eta );
			double dPcSF_a = eta * std::pow( 1.0 - a , -eta-1 );
			PcSF = PcSF_a + dPcSF_a * ( Sh - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::PcSF1( Sh ),   Sh = " << Sh << std::endl;
			exit(0);
		}
		return PcSF ;
	}

	double dPcSF1_dSh( double Sh ) const {
		double dPcSF=0. ;
		double eta = (lambda() * m - 1.)/(lambda()*m);
		double a = 0.95;
		if( Sh < a ){
			dPcSF = eta * std::pow( 1.0-Sh , -eta-1.0 );
		}
		else if( Sh >= a ){
			double dPcSF_a = eta * pow( 1.0 - a , -eta-1.0 );
			double ddPcSF_a = eta * (eta+1.) * std::pow( 1.0 - a , -eta-2.0 );
			dPcSF = dPcSF_a + ddPcSF_a * ( Sh - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::dPcSF1:dSh( Sh ),   Sh = " << Sh << std::endl;
			exit(0);
		}
		return dPcSF ;
	}

	double PcSF2( double phi_0, double phi ) const {
		double PcSF;
		double beta = betaf(); 
		double a = 0.05 ;
		if ( phi > a ){
			double term = ( phi_0/phi ) * ( ( 1-phi )/( 1. - phi_0 ));
			PcSF =pow( term, beta );
		}
		else if( phi <= a ){
			double term_a = ( phi_0/a ) * ( ( 1-a )/( 1. - phi_0 ));
			double PcSF_a = std::pow( term_a,beta );
			double dPcSF_a = 0.;
			double C = std::pow( phi_0/( 1. - phi_0 ) , beta );
			dPcSF_a -= beta * C * std::pow( 1.-a , beta-1. ) * std::pow( a , -beta-1. );
			PcSF = PcSF_a + dPcSF_a * ( phi - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::PcSF2( phi_0, phi )" << std::endl;
			exit(0);
		}

		return PcSF ;
	}

	double dPcSF2_dPhi( double phi_0, double phi ) const {

		double dPcSF = 0.;
		double a = 0.05 ;
		double beta = betaf(); 
		if ( phi > a ){
			double C = std::pow( phi_0/( 1. - phi_0 ) , beta );
			dPcSF -= beta * C * std::pow( 1.-phi , beta-1. ) * std::pow( phi , -beta-1. );
		}
		else if( phi <= a ){
			double C = std::pow( phi_0/( 1. - phi_0 ) , beta );
			double dPcSF_a = 0.;
			dPcSF_a -= beta * C * std::pow( 1.-a , beta-1. ) * std::pow( a , -beta-1. );
			double ddPcSF_a = 0.;
			ddPcSF_a = beta * C * std::pow( 1.-a , beta-1. ) * std::pow( a , -beta-1. );
			ddPcSF_a *= ( (beta+1.)/a + (beta-1.)/(1.-a) );
			dPcSF = dPcSF_a + ddPcSF_a * ( phi - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::PcSF2( phi_0, phi )" << std::endl;
			exit(0);
		}

		return dPcSF ;

	}

	double krW( double Sw, double Sh ) const {

		double Swe = effectiveSw(Sw,Sh);

		return std::pow(Swe, (2.0/lambda() + 3.0) ) ;
	}

	double dkrW_dSwe( double Sw, double Sh ) const {

		double Swe = effectiveSw(Sw,Sh);

		return (2.0/lambda() + 3.0) * std::pow(Swe, (2.0/lambda() + 3.0)-1. ) ;
	}

	double krNW( double Sw, double Sh ) const {

		double Swe = effectiveSw(Sw,Sh);

		return std::pow(1.0-Swe, 2.0) * ( 1.0 - std::pow(Swe, (2.0/lambda() + 1.0) ) );
	}

	double dkrNW_dSwe( double Sw, double Sh ) const {

		double Swe = effectiveSw(Sw,Sh);

		return (-1.) * ( 2.*(1.-Swe) + (2./lambda() + 1.) * std::pow( Swe, 2.0/lambda() ) );
	}


	/* PERMEABILITY SCALING FACTORS */

	double KSF1( double Sh ) const {

		return std::pow( (1.0 - Sh) , (5*m + 4)/(2*m) );
	}

	double dKSF1_dSh( double Sh ) const {
		return (-1.) * (5*m + 4)/(2*m) * std::pow( (1.0 - Sh) , (5*m + 4)/(2*m)-1. );
	}

	double KSF2( double phi_0, double phi ) const {
		// POWER LAW MODEL PROPOSED BY CIVAN (2001)
		// Read " kinetic simulation of methane hydrate formation and issociation in porous media " by Xuefei Sun and Kishore Mohanty
		// phi_0 = basePorosity_initial (and NOT sediemntPorosity!! )
		// phi is basePorosity at current time
		double beta = betaf(); 
		double term1 = phi / phi_0 ;
		double term2 = ( 1-phi_0 ) / ( 1. - phi ) ;

		double KSF=0.;
		double a = 0.95 ;

		if( phi < a ){
			KSF = term1 * std::pow( ( term1 * term2 ) , 2*beta );
		}
		else if ( phi >= a ){
			double C = std::pow( 1-phi_0 , 2*beta ) / std::pow( phi_0 , 2*beta +1 );
			double KSF_a = C * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta );
			double dKSF_a =    C * ( 2*beta +1 ) * std::pow( a, 2*beta   ) * std::pow( 1.-a , -2*beta      )
							 + C * ( -2*beta   ) * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta - 1. );
			KSF = KSF_a + dKSF_a * ( phi - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::KSF2( phi_0, phi )" << std::endl;
			exit(0);
		}
		return KSF;
	}

	double dKSF2_dPhi( double phi_0, double phi ) const {
		std::cout<< " ERROR in HydraulicProperties::dKSF2_dPhi( phi_0, phi )" << std::endl;
		exit(0);
	}


};




 #endif /* PARAMETERS_HYDRAULICPROPERTIES_HH_ */



// template<typename GV,typename PTree>
// class HydraulicProperties
// {
// private:
// 	  const GV& gv;
// 	  const PTree& ptree;

// 	  const static int dim = GV::dimension;

// 	  const static int numOfParams  = 5;
// 	  const static int id_Pentry 	= 0;
// 	  const static int id_lambda 	= 1;
// 	  const static int id_Swr 		= 2;
// 	  const static int id_Sgr 		= 3;
// 	  const static int id_beta		= 4;

// 	  Parameters<PTree> parameter;
// 	  MeshParameters<PTree> mesh;
// 	  //Soil<GV,PTree> soil;
// 	  CharacteristicValues characteristicValue;

// public:

//   //! construct from grid view
//   HydraulicProperties (const GV& gv_,const PTree& ptree_)
// : gv( gv_ ),
//   ptree(ptree_),
//   parameter(ptree_),
//   mesh(ptree_)
// {}

// 	/* PARAMETERS FOR THE HYDRAULIC PROPERTY CONSTITUTIVE LAW (BROOKS-COREY) */
// 	std::vector<double>
// 	BrooksCoreyParameters( const typename GV::Traits::template Codim<0>::Entity& element,
// 	   	   	 	 	 	   const Dune::FieldVector<double,dim>& xlocal ) const {

// 		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

// 		std::vector<double> BCParams (numOfParams,0.);

// 		auto prop_L = parameter.layer_properties();

// 		BCParams[id_Pentry] = prop_L[0][2] ; /*Pa*/
// 		BCParams[id_lambda] = prop_L[0][3] ;
// 		BCParams[id_Swr]    = prop_L[0][4] ;
// 		BCParams[id_Sgr]    = prop_L[0][5] ;
// 		BCParams[id_beta]   = prop_L[0][6] ;

// 		if( mesh.isLens(x) ){
// 			BCParams[id_Pentry] = prop_L[1][2] ; /*Pa*/
// 			BCParams[id_lambda] = prop_L[1][3] ;
// 			BCParams[id_Swr]    = prop_L[1][4] ;
// 			BCParams[id_Sgr]    = prop_L[1][5] ;
// 			BCParams[id_beta]   = prop_L[1][6] ;
// 		}

// 		return BCParams;
// 	}

// 	/* EFFECTIVE SATURATION */


// 	double dSwe_dSw( double Sw,
// 					 double Swr,
// 					 double Sgr) const {

// 		double dSwe =  1./( 1. - Sgr - Swr );

// 		return dSwe;
// 	}


	

// 	double dPc_dSwe( double Swe,
// 					 double Pentry, /*Pa*/
// 					 double lambda ) const {

// 		double eta = (1/lambda);
// 		double dPc = 0.; /*Pa*/
// 		double a = 0.05 ;

// 		if( Swe > a ){
// 			dPc/*Pa*/ = Pentry * (-1./lambda) * std::pow( Swe , -(1./lambda) - 1. );
// 		}
// 		else if ( Swe <= a ){
// 			double dPc_a  = Pentry * (-1./lambda) * std::pow( a , -(1./lambda) - 1. ) ;
// 			double ddPc_a = Pentry * (-1./lambda) * (-1./lambda-1.) * std::pow( a , -(1./lambda) - 2. );
// 			dPc/*Pa*/ = dPc_a + ddPc_a * ( Swe - a );
// 		}
// 		else {
// 			std::cout<< " ERROR in HydraulicProperties::dPc_dSwe( Swe, Pentry, lambda ) "
// 					 << "  , Swe = " << Swe
// 					 << "  , dPc  = " << dPc << std::endl;
// //			exit(0);
// 		}

// 		return dPc; /*Pa*/
// 	}

// 	/* Pc SCALING FACTORS */

// 	double PcSF( double phi,
// 				 double phi_0,
// 				 double beta) const {

// 		double pcsf = 0.;

// 		double a = 0.05 ;
// 		if ( phi > a ){
// 			double term = ( phi_0/phi ) * ( ( 1-phi )/( 1. - phi_0 ));
// 			pcsf =pow( term, beta );
// 		}
// 		else if( phi <= a ){
// 			double term_a = ( phi_0/a ) * ( ( 1-a )/( 1. - phi_0 ));
// 			double pcsf_a = std::pow( term_a,beta );
// 			double dpcsf_a = 0.;
// 			double C = std::pow( phi_0/( 1. - phi_0 ) , beta );
// 			dpcsf_a -= beta * C * std::pow( 1.-a , beta-1. ) * std::pow( a , -beta-1. );
// 			pcsf = pcsf_a + dpcsf_a * ( phi - a );
// 		}
// 		else {
// 			std::cout<< " ERROR in " << __FILE__
// 					 << " function: PcSF( phi, phi_0, beta )" << std::endl;
// //			exit(0);
// 		}

// 		return pcsf ;
// 	}

// 	/* RELATIVE PERMEABILITIES */

// 	double krW( const typename GV::Traits::template Codim<0>::Entity& element,
// 	   	 	 	const Dune::FieldVector<double,dim>& xlocal ,
// 				double Sw ) const {

// 		auto BCParams = BrooksCoreyParameters(element,xlocal);
// 		double lambda 	= BCParams[id_lambda];
// 		double Sgr 		= BCParams[id_Sgr];
// 		double Swr 		= BCParams[id_Swr];

// 		double Swe = EffectiveSw(Sw,Swr,Sgr);

// 		double kr = std::pow(Swe, (2.0/lambda + 3.0) );
// 		if( Swe>1.){
// 			kr=1.;
// 		}
// 		if( Swe<0.){
// 			kr=0.;
// 		}

// 		return kr ;
// 	}

// 	double krN( const typename GV::Traits::template Codim<0>::Entity& element,
// 	   	 	 	const Dune::FieldVector<double,dim>& xlocal ,
// 				double Sw ) const {

// 		auto BCParams = BrooksCoreyParameters(element,xlocal);
// 		double lambda 	= BCParams[id_lambda];
// 		double Sgr 		= BCParams[id_Sgr];
// 		double Swr 		= BCParams[id_Swr];

// 		double Swe = EffectiveSw(Sw,Swr,Sgr);

// 		double kr = std::pow(1.0-Swe, 2.0) * ( 1.0 - std::pow(Swe, (2.0/lambda + 1.0) ) );

// 		if( Swe>1.){
// 			kr=0.;
// 		}
// 		if( Swe<0.){
// 			kr=1.;
// 		}
// 		return kr;
// 	}

// 	/* PERMEABILITY SCALING FACTORS */

// 	double PermeabilityScalingFactor( const typename GV::Traits::template Codim<0>::Entity& element,
// 	   	 	 	 	 	 	 	 	  const Dune::FieldVector<double,dim>& xlocal ,
// 									  double porosity ) const {

// 		auto BCParams = BrooksCoreyParameters(element,xlocal);
// 		double Pentry 	= BCParams[id_Pentry];
// 		double lambda 	= BCParams[id_lambda];
// 		double Sgr 		= BCParams[id_Sgr];
// 		double Swr 		= BCParams[id_Swr];
// 		double beta 	= BCParams[id_beta];

// 		double porosity_0 = 0.3;//soil.SedimentPorosity( element,xlocal );
// 		double SF_por = KSF( porosity, porosity_0, beta );

// 		double SF =SF_por;
// 		return SF;
// 	}

// 	double KSF( double phi,
// 				double phi_0,
// 				double beta ) const {
// 		// POWER LAW MODEL PROPOSED BY CIVAN (2001)
// 		// Read " kinetic simulation of methane hydrate formation and issociation in porous media " by Xuefei Sun and Kishore Mohanty
// 		// phi_0 = basePorosity_initial (and NOT sediemntPorosity!! )
// 		// phi is basePorosity at current time

// 		double term1 = phi / phi_0 ;
// 		double term2 = ( 1-phi_0 ) / ( 1. - phi ) ;

// 		double ksf=0.;
// 		double a = 0.95 ;

// 		if( phi < a ){
// 			ksf = term1 * std::pow( ( term1 * term2 ) , 2*beta );
// 		}
// 		else if ( phi >= a ){
// 			double C = std::pow( 1-phi_0 , 2*beta ) / std::pow( phi_0 , 2*beta +1 );
// 			double ksf_a = C * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta );
// 			double dksf_a =    C * ( 2*beta +1 ) * std::pow( a, 2*beta   ) * std::pow( 1.-a , -2*beta      )
// 							 + C * ( -2*beta   ) * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta - 1. );
// 			ksf = ksf_a + dksf_a * ( phi - a );
// 		}
// 		else {
// 			std::cout<< " ERROR in " << __FILE__
// 					 << " function: KSF2( phi, phi_0, beta )" << std::endl;
// //			exit(0);
// 		}
// 		return ksf;
// 	}


//   //! get a reference to the grid view
//   inline const GV& getGridView () {return gv;}

// };

/*
 * problem_20160922_test1.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PROBLEMSPECS_PROBLEM_20160922_TEST1_HH_
#define PROBLEMSPECS_PROBLEM_20160922_TEST1_HH_

#include"../parameters/hydraulicProperties.hh"

/* PROBLEM DEFINITION:
 *
 *
 */

class Problem_20160922_test1{
private:
	constexpr static double pi = atan(1.)*4 ;
	constexpr static double eps = 1.0e-6;
	Indices indices;
	CharacteristicValues characteristicValues;

	const static int refLevel = 1;
public:

	/******************************************************************************/

	std::string getPathName(){

		std::string pathName = "/home/amir/dune-master/proj-HydrateDG/dune/proj-HydrateDG/hydrate-DG-old/outputs/";
	    return pathName ;
	}

	std::string getFileName(){

		std::string test = "test_";
	    return test;

	}

	/*
	 * 2D -> X and Z
	 * 3D -> X, Y and Z
	 */
	const static int dimension 		= 2  	;
	constexpr static double origin 		= 0.0	;
	constexpr static double X_length 	= 1	;
	constexpr static double Z_length  	= 1.0	;
	constexpr static double Y_length	= 1.	;		//only if dim=3

	const static int X_cells	= 8 ;
	const static int Z_cells	= 8 ;
	const static int Y_cells	= 1	 ;
	/**********************************************************************/

	constexpr static double T = 283/CharacteristicValues::T_c; /*K*/


	double SedimentPorosity( Dune::FieldVector< double, dimension > globalPos ) const {
		double porosity = 0.4;

		return porosity ;

	}

	double SedimentPermeability( Dune::FieldVector< double, dimension > globalPos ) const {

		double permeability = 1. * 1.e-14;
		return permeability ;

	}

	Dune::FieldMatrix<double,dimension,dimension> SedimentPermeabilityTensor( Dune::FieldVector< double, dimension > globalPos ) const {

		double K = SedimentPermeability(globalPos);
		Dune::FieldMatrix<double,dimension,dimension> PermeabilityTensor;

		for (std::size_t i=0; i<dimension; i++){
			for (std::size_t j=0; j<dimension; j++){
				PermeabilityTensor[i][j] = (i==j) ? K : 0;
			}
		}

		return PermeabilityTensor;

	}

	double SoilGrainRadius( Dune::FieldVector< double, dimension > globalPos ) const {
	  // Bear et at, 1972
	  double rp = sqrt( 45.0 * SedimentPermeability( globalPos ) * pow( 1- SedimentPorosity( globalPos ) , 2.0 )/pow( SedimentPorosity( globalPos ),3.0) );
// 	  std::cout<< "rp = " << rp << std::endl;
	      return rp;
	}


	double AreaFactor() const {
	   return 1.;
	}

	double ReactionRateParam_k0() const {

		double k0 = 360000.;		//defined in mol/m².Pa.s

		return k0;
	}

	double ReactionRateParam_deltaEabyR() const {
		return 9752.73;
	}

	double ReactionRateParam_kForm() const {
		//Englezos et al. (1987a)			reference article: Kinetic simulation of methane hydrate formation and dissociation in porous media
											//  ------------------   author: K Mohanty
		double kd_form ;		//defined in mol/m².Pa.s
		kd_form = 0.;//0.5875 * 1.e-11;

		return kd_form;
	}

	/******************************************************************************/

	bool isLeftBoundary( Dune::FieldVector< double, dimension > globalPos ) const{
		if( globalPos[0] < origin + eps )
			return true;
		else
			return false;
	}

	bool isRightBoundary( Dune::FieldVector< double, dimension > globalPos ) const{
		if( globalPos[0] > X_length - eps )
			return true;
		else
			return false;
	}

	bool isBottomBoundary( Dune::FieldVector< double, dimension > globalPos ) const{
		if( dimension == 2 ){
			if( globalPos[1] < origin + eps ){
				return true;
			}
			else
				return false;
		}
		else if( dimension == 3 ){
			if( globalPos[2] < origin + eps ){
				return true;
			}
			else
				return false;
		}
	}

	bool isTopBoundary( Dune::FieldVector< double, dimension > globalPos ) const{
		if( dimension == 2 ){
			if( globalPos[1] > Z_length - eps ){
				return true;
			}
			else
				return false;
		}
		else if( dimension == 3 ){
			if( globalPos[2] < Z_length - eps ){
				return true;
			}
			else
				return false;
		}
	}

	bool isFrontBoundary( Dune::FieldVector< double, dimension > globalPos ) const{
		if( dimension == 3 ){
			if( globalPos[1] < origin + eps ){
				return true;
			}
			else
				return false;
		}
		else
			return false;
	}

	bool isBackBoundary( Dune::FieldVector< double, dimension > globalPos ) const{
		if( dimension == 3 ){
			if( globalPos[1] > Y_length - eps ){
				return true;
			}
			else
				return false;
		}
		else
			return false;
	}

	/**********************************************************************/
	std::vector< double > ProblemICValues ( Dune::FieldVector< double, dimension > globalPos ) const {

		std::vector< double > icvalue(Indices::numOfFlowPVs);
		for(int i = 0; i<Indices::numOfFlowPVs ; i++){
			icvalue[i] = 0.; //
		}

		double Pg = 10.*1.e6; /*Pa*/
		double Sw = 0.3;
		double Sh = 0.4;
		double T = 273.15+10.; /*K*/
		double XCH4 = 0.1;
		double YH2O = 0.1;
		double XC = 5.5e-3;
		HydraulicProperties hydraulicProperty;
		double Pc = hydraulicProperty.suctionPressure(Sw,Sh) * hydraulicProperty.PcSF1(Sh);

		icvalue[Indices::PVId_Pg] = Pg ; //P_w + P_c ;
		icvalue[Indices::PVId_Sw] = Sw ;
		icvalue[Indices::PVId_Sh] = Sh ;
		icvalue[Indices::PVId_Pc] = Pc ;
		icvalue[Indices::PVId_T ] = T  ;
		icvalue[Indices::PVId_XCH4] = XCH4 ;
		icvalue[Indices::PVId_YH2O ] = YH2O  ;
		icvalue[Indices::PVId_C] = XC ;

		return icvalue;
	}

	/******************************************************************************/

	std::vector< int > ProblemBCTypes ( Dune::FieldVector< double, dimension > globalPos ) const
	{

		std::vector< int > bctype(Indices::numOfFlowPVs);
		for(int i = 0; i<Indices::numOfFlowPVs ; i++){
			bctype[i] = 0; //
		}

		if( isTopBoundary(globalPos) ){
			bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			bctype[indices.PVId_T ] = indices.BCId_neumann ;
			bctype[indices.PVId_XCH4] = indices.BCId_neumann ;
			bctype[indices.PVId_YH2O ] = indices.BCId_neumann ;
			bctype[indices.PVId_C] = indices.BCId_neumann ;
		}
		else if( isBottomBoundary(globalPos) ){

			bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			bctype[indices.PVId_T ] = indices.BCId_neumann ;
			bctype[indices.PVId_XCH4] = indices.BCId_neumann ;
			bctype[indices.PVId_YH2O ] = indices.BCId_neumann ;
		}
		else if( isLeftBoundary(globalPos) ){
			//if( (globalPos[1] > Z_length*(3./8.) - eps) and (globalPos[1] < Z_length*(5./8.) + eps) ){
				bctype[indices.PVId_Pg] = indices.BCId_dirichlet ;
				bctype[indices.PVId_Pc] = indices.BCId_dirichlet ;
				bctype[indices.PVId_Sw] = indices.BCId_dirichlet ;
			// }else{
			// 	bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			// 	bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			// 	bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			// }
			bctype[indices.PVId_T ] = indices.BCId_neumann ;
			bctype[indices.PVId_XCH4] = indices.BCId_dirichlet ;
			bctype[indices.PVId_YH2O ] = indices.BCId_dirichlet ;
			bctype[indices.PVId_C] = indices.BCId_neumann ;
		}
		else if( isRightBoundary(globalPos) ){
			bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			bctype[indices.PVId_T ] = indices.BCId_neumann ;
			bctype[indices.PVId_XCH4] = indices.BCId_neumann ;
			bctype[indices.PVId_YH2O ] = indices.BCId_neumann ;
			bctype[indices.PVId_C] = indices.BCId_neumann ;
		}

		return bctype;
	}

	/******************************************************************************/

	std::vector< double > ProblemBCValues ( Dune::FieldVector< double, dimension > globalPos , double time ) const {

		std::vector< double > bcvalue(Indices::numOfFlowPVs);
		for(int i = 0; i<Indices::numOfFlowPVs ; i++){
			bcvalue[i] = 0.; //
		}

		if( isTopBoundary(globalPos) ){
			bcvalue[indices.PVId_Pg] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Sw] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_T ] = 0.;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
			bcvalue[indices.PVId_XCH4] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_YH2O ] = 0.;
			bcvalue[indices.PVId_C] = 0.;
		}
		else if( isBottomBoundary(globalPos) ){
			bcvalue[indices.PVId_Pg] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Sw] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_T ] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_XCH4] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_YH2O ] = 0.;
		}
		else if( isLeftBoundary(globalPos) ){
			HydraulicProperties hydraulicProperty;
			//if( (globalPos[1] > Z_length*(3./8.) - eps) and (globalPos[1] < Z_length*(5./8.) + eps) ){
				bcvalue[indices.PVId_Pg] = 6.*1.e6;//11.*1.e6;//indices.BCId_dirichlet ;
				bcvalue[indices.PVId_Pc] = hydraulicProperty.Pentry();//indices.BCId_neumann ;
				bcvalue[indices.PVId_Sw] = 0.;//indices.BCId_neumann ;
			//}else{
			// 	bcvalue[indices.PVId_Pg] = 0.;//indices.BCId_neumann ;
			// 	bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
			// 	bcvalue[indices.PVId_Sw] = 0.;//indices.BCId_neumann ;
			// }
			bcvalue[indices.PVId_T ] = 0.1;//indices.BCId_neumann ;
			bcvalue[indices.PVId_XCH4] = 0.1;//indices.BCId_neumann ;
			bcvalue[indices.PVId_YH2O ] = 0.1;
			bcvalue[indices.PVId_C] = 0.;
		}
		else if( isRightBoundary(globalPos) ){
			bcvalue[indices.PVId_Pg] = 0.;//ProblemICValues(globalPos)[indices.PVId_Pg];//indices.BCId_dirichlet ;
			bcvalue[indices.PVId_Pc] = 0.;//ProblemICValues(globalPos)[indices.PVId_Pc];//indices.BCId_dirichlet ;
			bcvalue[indices.PVId_Sw] = 0.;//ProblemICValues(globalPos)[indices.PVId_Sw];//indices.BCId_dirichlet ;
			bcvalue[indices.PVId_T ] = 0.;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
			bcvalue[indices.PVId_XCH4] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_YH2O ] = 0.;
			bcvalue[indices.PVId_C] = 0.;
		}

		return bcvalue;
	}



	/******************************************************************************/

};

#endif /* PROBLEMSPECS_PROBLEM_20160922_TEST1_HH_ */

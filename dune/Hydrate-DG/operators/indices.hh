/*
 * indices.hh
 *
 * 			DONE
 */

class Indices{
private:

public:

	/* PHASES */
#ifdef PLOT_VELOCITIES
	static const unsigned int numOfPhases = 2;
	static const unsigned int phaseId_g = 0;
	static const unsigned int phaseId_w = 1;
#endif
	//static const unsigned int phaseId_liquid = 0;
	//static const unsigned int phaseId_gas 	 = 1;

	/* COMPONENTS */
	static const int numOfComps  = 4;
	static const int compId_XCH4 = 0;
	static const int compId_XH2O = 1;
	static const int compId_YCH4 = 2;
	static const int compId_YH2O = 3;


	/* PRIMARY VARIABLES */  // NOTE: The order of the indices must be the same as the order in the GFS, initial and boundary conditions
	static const int numOfPVs 	 = 7 ;
	static const int PVId_Pw  	= 0;
	static const int PVId_Sg  	= 1;
	static const int PVId_Sh	= 2;
	static const int PVId_T 	= 3;
	static const int PVId_XCH4	= 4 ;
	static const int PVId_YH2O	= 5 ;
	static const int PVId_C   	= 6;

	static const int numOfVs 	 = 5 ;
	static const int VId_Pw  	= 0;
	static const int VId_Sg  	= 1;
	static const int VId_XCH4	= 2 ;
	static const int VId_YH2O	= 3 ;
	static const int VId_XC   	= 4;

	static const unsigned int numOfSVs 	= 16;
	static const unsigned int SVId_Pg	= 0;  // gas phase pressure
	static const unsigned int SVId_Pw	= 1;  // water phase pressure
	static const unsigned int SVId_Pc	= 2;  // capillary pressure
	static const unsigned int SVId_Sg	= 3;  // gas saturation
	static const unsigned int SVId_Sw	= 4;  // water saturation
	static const unsigned int SVId_Sh	= 5;  // hydrate saturation
	static const unsigned int SVId_XCH4	= 6;  // CH4 mole fraction in water
	static const unsigned int SVId_XH2O	= 7;  // H2O mole fraction in water
	static const unsigned int SVId_YCH4	= 8;  // CH4 mole fraction in gas
	static const unsigned int SVId_YH2O	= 9;  // H2O mole fraction in gas
	static const unsigned int SVId_Xc	= 10; // salt mole fraction in water
	static const unsigned int SVId_T	= 11; // temperature
 	static const unsigned int SVId_K 	= 12; // absolute permeability
	static const unsigned int SVId_z	= 13; // methane gas compressibility factor
	static const unsigned int SVId_Peq 	= 14; // hydrate equilibrium pressure
	static const unsigned int SVId_por	= 15; // total porosity

	/* SECONDARY VARIABLES FOR POST PROCESS AND O/P */
	static const unsigned int numOfFlowPPvars = 11 ;
	static const unsigned int PP_Sw		= 0 ;
	static const unsigned int PP_Pg		= 1 ;
	static const unsigned int PP_XCH4	= 2 ;
	static const unsigned int PP_XH2O	= 3 ;
	static const unsigned int PP_YCH4	= 4 ;
	static const unsigned int PP_YH2O	= 5 ;
	static const unsigned int PP_Pc		= 6 ;
	static const unsigned int PP_K		= 7 ;
	static const unsigned int PP_krw	= 8 ;
	static const unsigned int PP_krg	= 9 ;
	static const unsigned int PP_zCH4	= 10;

	/* BOUNDARY CONDITION */
	static const unsigned int numOfVelBCs 	= 4; // V_g , V_w
	static const unsigned int BCId_water	= 0;
	static const unsigned int BCId_gas		= 1;
	static const unsigned int BCId_salt 	= 2;
	static const unsigned int BCId_heat		= 3;
	static const unsigned int BCId_dirichlet		= 1;
	static const unsigned int BCId_neumann 	 		= 0;
	static const unsigned int BCId_depressurization = 2;



};



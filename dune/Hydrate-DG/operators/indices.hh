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
	// static const int PVId_pen   = 7;

	static const int numOfVs 	 = 6 ;
	static const int VId_Pw  	= 0;
	static const int VId_Sg  	= 1;
	static const int VId_Sh		= 2;
	static const int VId_T 		= 3;
	static const int VId_XCH4	= 4 ;
	static const int VId_YH2O	= 5 ;

	const static int numOfSVs 	= 32;
	const static int SVId_Pg	= 0;  // gas phase pressure
	const static int SVId_Pw	= 1;  // water phase pressure
	const static int SVId_Pc	= 2;  // capillary pressure
	const static int SVId_Sg	= 3;  // gas saturation
	const static int SVId_Sh	= 4;  // hydrate saturation
	const static int SVId_Sw	= 5;  // water saturation
	const static int SVId_T		= 6;  // Temp
	const static int SVId_XCH4	= 7;  // CH4 mole fraction in water
	const static int SVId_XC	= 8;  // CH4 mole fraction in water
	const static int SVId_XH2O	= 9;  // H2O mole fraction in water
	const static int SVId_YCH4	= 10;  // CH4 mole fraction in gas
	const static int SVId_YH2O	= 11;  // H2O mole fraction in gas
 	const static int SVId_rhow 	= 12; // water density
 	const static int SVId_rhog 	= 13; // methane gas density
 	const static int SVId_K 	= 14; // absolute permeability
 	const static int SVId_krw 	= 15; // rel. water permeability
 	const static int SVId_krg 	= 16; // rel. gas permeability
 	const static int SVId_muw 	= 17; // dyn. water viscosity
 	const static int SVId_mug 	= 18; // dyn. gas viscosity
	const static int SVId_zCH4	= 19; // methane gas compressibility factor
	const static int SVId_por	= 20; // total porosity
 	const static int SVId_DH2O 	= 21; // H2O binary diffusion coeff. in methane
 	const static int SVId_DCH4 	= 22; // CH4 binary diffusion coeff. in water
 	const static int SVId_Pwsat = 23; // Saturation pressure for water vapor
 	const static int SVId_HCH4 	= 24; // Henry's constant
 	const static int SVId_tau	= 25; // tortuosity
 	const static int SVId_Peq	= 26; // Equi. Pressure
 	const static int SVId_HS	= 27; // Hydrate Stability zone
 	const static int SVId_Vwx	= 28; // Velocity
	const static int SVId_Vwy	= 29; // Velocity
	const static int SVId_Vgx	= 30; // Velocity
	const static int SVId_Vgy	= 31; // Velocity


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



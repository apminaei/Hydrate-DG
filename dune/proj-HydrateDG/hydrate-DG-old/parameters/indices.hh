/*
 * indices.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_INDICES_HH_
#define PARAMETERS_INDICES_HH_

class Indices{
private:

public:

	/* PHASES */
	static const unsigned int phaseId_liquid = 0;
	static const unsigned int phaseId_gas 	 = 1;

	/* COMPONENTS */
	static const unsigned int compId_H2O 	 = 0;
	static const unsigned int compId_CH4 	 = 1;


	/* PRIMARY VARIABLES */
	static const unsigned int numOfFlowPVs 	 = 8 ;
	static const unsigned int PVId_Pg  	= 0;
	static const unsigned int PVId_Pc  	= 1;
	static const unsigned int PVId_Sw  	= 2;
	static const unsigned int PVId_Sh 	= 3;
	static const unsigned int PVId_T   	= 4;
	static const unsigned int PVId_XCH4	= 5 ;
	static const unsigned int PVId_YH2O	= 6 ;
	static const unsigned int PVId_C   	= 7;

	/* SECONDARY VARIABLES FOR POST PROCESS AND O/P */
	static const unsigned int numOfFlowPPvars = 11 ;
	static const unsigned int PP_Sg		= 0 ;
	static const unsigned int PP_Pw		= 1 ;
	static const unsigned int PP_XCH4	= 2 ;
	static const unsigned int PP_XH2O	= 3 ;
	static const unsigned int PP_YCH4	= 4 ;
	static const unsigned int PP_YH2O	= 5 ;
	static const unsigned int PP_Pc		= 6 ;
	static const unsigned int PP_K		= 7 ;
	static const unsigned int PP_krw	= 8 ;
	static const unsigned int PP_krg	= 9 ;
	static const unsigned int PP_zCH4	= 10;

	/* BOUNDARY CONDITION TYPES */
	static const unsigned int BCId_dirichlet 		= 0;
	static const unsigned int BCId_neumann 	 		= 1;
	static const unsigned int BCId_depressurization = 2;

};



#endif /* PARAMETERS_INDICES_HH_ */
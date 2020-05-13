template<typename GV,typename Properties>
class ProblemBoundaryConditions
{
private :
	const GV& gv ;
	const Properties& property;
	const static int dim = GV::dimension;
	Indices indices;
	CharacteristicValues characteristicValues;

public :

	// ! construct from gridview
	ProblemBoundaryConditions (const GV& gv_,const Properties& property_)
	: gv ( gv_ ),
	  property(property_)
	{}
	
	/* NOTE:
	 * dirichlet: BC_w -> Pw, BC_g -> Sg, 
	 * neumann: total mass flux (convective+diffusive) BC_w -> f_w, BC_g -> f_g
	 * */
	
	/* boundary types */
	template<typename I> std::vector< int >
	type( I& intersection,
		  const Dune::FieldVector<double,dim>& globalPos,
		  double time/*s*/,
		  double dt/*s*/ ) const {
		
		//auto xglobal = intersection.geometry().global(xlocal);
		
		std::vector< int > bctype(Indices::numOfPVs, 0);

		bctype[indices.PVId_T ] = indices.BCId_dirichlet ;
		bctype[indices.PVId_C ] = indices.BCId_neumann ;
		bctype[indices.PVId_Sh] = indices.BCId_neumann ;
		//bctype[indices.PVId_Pc] = indices.BCId_neumann ;
		bctype[indices.PVId_Sg] = indices.BCId_neumann ;
		bctype[indices.PVId_XCH4 ] = indices.BCId_neumann ;
		bctype[indices.PVId_YH2O ] = indices.BCId_neumann ;
		//if( property.mesh.isTopBoundary(globalPos) ){
		if (0. < time <= 200./* hrs */) {
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		}
		else if ( 200. < time <= 350./* hrs */) {
			bctype[indices.PVId_Pw] = indices.BCId_neumann ;
		}
		else {
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		}
		
		return bctype;
	}

	/* boundary values */
	template<typename I> std::vector< double >
	value ( I& intersection,
			const Dune::FieldVector<double,dim>& globalPos,
			double time/*s*/,
			double dt/*s*/ ) const {

		//auto xglobal = intersection.geometry().global(xlocal);

		std::vector< double > bcvalue(Indices::numOfPVs,0.);
		
		bcvalue[indices.PVId_T ] = 4.+273.15 ;
		bcvalue[indices.PVId_C ] = 0. ;
		bcvalue[indices.PVId_Sh] = 0. ;
		//bcvalue[indices.PVId_Pc] = 0. ;
		bcvalue[indices.PVId_Sg] = 0. ;
		bcvalue[indices.PVId_XCH4 ] = 0. ;
		bcvalue[indices.PVId_YH2O ] = 0. ;
		//if( property.mesh.isTopBoundary(globalPos) ){
		if (0. < time <= 200./* hrs */) {
			bcvalue[indices.PVId_Pw] = 2.e6 ; /* Pa */
		}
		else if ( 200. < time <= 350./* hrs */) {
			bcvalue[indices.PVId_Pw] = 0. ;
		}
		else if ( 350. < time /* hrs */) {
			bcvalue[indices.PVId_Pw] = 5.e6 ;
		}
		else{
			bcvalue[indices.PVId_Pw] = 5.e6 + 10. * (time - 450 * 3600) ;
		}

		// if( property.mesh.isTopBoundary(globalPos) ){
		// 	bcvalue[indices.PVId_Pw] = 0.;//indices.BCId_neumann ;
		// 	bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
		// 	bcvalue[indices.PVId_Sg] = 0.;//indices.BCId_neumann ;
		// 	//bcvalue[indices.PVId_T ] = 0.;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
		// }
		// else if( property.mesh.isBottomBoundary(globalPos) ){
		// 	bcvalue[indices.PVId_Pw] = 0.;//indices.BCId_neumann ;
		// 	bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
		// 	bcvalue[indices.PVId_Sg] = 0.;//indices.BCId_neumann ;
		// 	//bcvalue[indices.PVId_T ] = 0.;//indices.BCId_neumann ;
		// }
		// else if( property.mesh.isLeftBoundary(globalPos) ){
		// 	HydraulicProperties hydraulicProperty;
		// 	//if( (globalPos[1] > Z_length*(3./8.) - eps) and (globalPos[1] < Z_length*(5./8.) + eps) ){
		// 		bcvalue[indices.PVId_Pw] = 6.*1.e6;//11.*1.e6;//indices.BCId_dirichlet ;
		// 		bcvalue[indices.PVId_Pc] = hydraulicProperty.Pentry();//indices.BCId_neumann ;
		// 		bcvalue[indices.PVId_Sg] = 0.;//indices.BCId_neumann ;
		// 	//}else{
		// 	// 	bcvalue[indices.PVId_Pw] = 0.;//indices.BCId_neumann ;
		// 	// 	bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
		// 	// 	bcvalue[indices.PVId_Sg] = 0.;//indices.BCId_neumann ;
		// 	// }
		// 	//bcvalue[indices.PVId_T ] = 0.;//indices.BCId_neumann ;
		// }
		// else if( property.mesh.isRightBoundary(globalPos) ){
		// 	bcvalue[indices.PVId_Pw] = 0.;//ProblemICValues(globalPos)[indices.PVId_Pw];//indices.BCId_dirichlet ;
		// 	bcvalue[indices.PVId_Pc] = 0.;//ProblemICValues(globalPos)[indices.PVId_Pc];//indices.BCId_dirichlet ;
		// 	bcvalue[indices.PVId_Sg] = 0.;//ProblemICValues(globalPos)[indices.PVId_Sg];//indices.BCId_dirichlet ;
		// 	//bcvalue[indices.PVId_T ] = 0.;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
		// }


		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};

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
		if( property.mesh.isTopBoundary(globalPos) ){
			bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			//bctype[indices.PVId_T ] = indices.BCId_neumann ;
		}
		else if( property.mesh.isBottomBoundary(globalPos) ){

			bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			//bctype[indices.PVId_T ] = indices.BCId_neumann ;
		}
		else if( property.mesh.isLeftBoundary(globalPos) ){
			//if( (globalPos[1] > Z_length*(3./8.) - eps) and (globalPos[1] < Z_length*(5./8.) + eps) ){
				bctype[indices.PVId_Pg] = indices.BCId_dirichlet ;
				bctype[indices.PVId_Pc] = indices.BCId_dirichlet ;
				bctype[indices.PVId_Sw] = indices.BCId_dirichlet ;
			// }else{
			// 	bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			// 	bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			// 	bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			// }
			//bctype[indices.PVId_T ] = indices.BCId_neumann ;
		}
		else if( property.mesh.isRightBoundary(globalPos) ){
			bctype[indices.PVId_Pg] = indices.BCId_neumann ;
			bctype[indices.PVId_Pc] = indices.BCId_neumann ;
			bctype[indices.PVId_Sw] = indices.BCId_neumann ;
			//bctype[indices.PVId_T ] = indices.BCId_neumann ;
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
		// if( property.mesh.isLeftBoundary(xglobal) ){
		// 	bcv[indices.PVId_Pg] = property.parameter.LeftPw(xglobal);
		// 	if( property.mesh.isInlet(xglobal) ) {
		// 		bcv[indices.PVId_Pg] = property.parameter.InletSg(xglobal) ;
		// 	}else bcv[indices.PVId_Pg] = property.parameter.InitialSg(xglobal) ;
		// }

		// if( property.mesh.isRightBoundary(xglobal) ){
		// 	bcv[indices.PVId_Pg] = property.parameter.InitialPw(xglobal);
		// 	bcv[indices.PVId_Pg] = property.parameter.InitialSg(xglobal) ;
		// }

		if( property.mesh.isTopBoundary(globalPos) ){
			bcvalue[indices.PVId_Pg] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Sw] = 0.;//indices.BCId_neumann ;
			//bcvalue[indices.PVId_T ] = 0.;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
		}
		else if( property.mesh.isBottomBoundary(globalPos) ){
			bcvalue[indices.PVId_Pg] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Pc] = 0.;//indices.BCId_neumann ;
			bcvalue[indices.PVId_Sw] = 0.;//indices.BCId_neumann ;
			//bcvalue[indices.PVId_T ] = 0.;//indices.BCId_neumann ;
		}
		else if( property.mesh.isLeftBoundary(globalPos) ){
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
			//bcvalue[indices.PVId_T ] = 0.;//indices.BCId_neumann ;
		}
		else if( property.mesh.isRightBoundary(globalPos) ){
			bcvalue[indices.PVId_Pg] = 0.;//ProblemICValues(globalPos)[indices.PVId_Pg];//indices.BCId_dirichlet ;
			bcvalue[indices.PVId_Pc] = 0.;//ProblemICValues(globalPos)[indices.PVId_Pc];//indices.BCId_dirichlet ;
			bcvalue[indices.PVId_Sw] = 0.;//ProblemICValues(globalPos)[indices.PVId_Sw];//indices.BCId_dirichlet ;
			//bcvalue[indices.PVId_T ] = 0.;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
		}


		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};

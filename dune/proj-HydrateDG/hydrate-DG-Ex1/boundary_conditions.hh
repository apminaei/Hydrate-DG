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
		  const Dune::FieldVector<double,dim>& xglobal,
		  double time/*s*/,
		  double dt/*s*/ ) const {
		
		//auto xglobal = intersection.geometry().global(xlocal);
		
		std::vector< int > bctype(Indices::numOfPVs, 0);
		if( property.mesh.isLeftBoundary(xglobal) or property.mesh.isRightBoundary(xglobal) ){
			bctype[indices.PVId_Pg] = indices.BCId_dirichlet ;
			bctype[indices.PVId_Pc] = indices.BCId_dirichlet ;
			bctype[indices.PVId_Sw] = indices.BCId_dirichlet ;
			
			bctype[indices.PVId_XCH4] = indices.BCId_dirichlet ;
			bctype[indices.PVId_YH2O ] = indices.BCId_dirichlet ;
			bctype[indices.PVId_C] = indices.BCId_neumann ;
		}else{
			bctype[indices.PVId_Pg] = indices.BCId_dirichlet ;
			bctype[indices.PVId_Pc] = indices.BCId_dirichlet ;
			bctype[indices.PVId_Sw] = indices.BCId_dirichlet ;
			
			bctype[indices.PVId_XCH4] = indices.BCId_dirichlet ;
			bctype[indices.PVId_YH2O ] = indices.BCId_dirichlet ;
			bctype[indices.PVId_C] = indices.BCId_neumann ;
		}
		
		return bctype;
	}

	/* boundary values */
	template<typename I> std::vector< double >
	value ( I& intersection,
			const Dune::FieldVector<double,dim>& xglobal,
			double time/*s*/,
			double dt/*s*/ ) const {

		//auto xglobal = intersection.geometry().global(xlocal);

		std::vector< double > bcv(Indices::numOfPVs,0.);
		if( property.mesh.isLeftBoundary(xglobal) ){
			bcv[indices.PVId_Pg] = property.parameter.LeftPw(xglobal);
			if( property.mesh.isInlet(xglobal) ) {
				bcv[indices.PVId_Pg] = property.parameter.InletSg(xglobal) ;
			}else bcv[indices.PVId_Pg] = property.parameter.InitialSg(xglobal) ;
		}

		if( property.mesh.isRightBoundary(xglobal) ){
			bcv[indices.PVId_Pg] = property.parameter.InitialPw(xglobal);
			bcv[indices.PVId_Pg] = property.parameter.InitialSg(xglobal) ;
		}

		return bcv;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};

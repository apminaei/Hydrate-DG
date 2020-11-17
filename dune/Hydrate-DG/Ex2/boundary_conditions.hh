/* ALL VALUES ARE NONDIMENSIONAL */
template<typename GV,typename Properties>
class ProblemBoundaryConditions
{
private :
	const GV& gv ;
	const Properties& property;
	const static int dim = GV::dimension;
	Indices indices;
	CharacteristicValues characteristicValues;
	ProblemInitialConditions<GV,Properties> icvalue;
	// double time_fraction = property.parameter.time_end() / 31.536e6; 
	double Xc_time = 1. / (364*24*3600);
	double press_rate = 1.;

public :

	// ! construct from gridview
	ProblemBoundaryConditions (const GV& gv_,const Properties& property_)
	: gv ( gv_ ),
	  property(property_),
	  icvalue(gv_, property_)
	{}
	
	/* NOTE:
	 * dirichlet: BC_w -> Pw, BC_g -> Sg, 
	 * neumann: total mass flux (convective+diffusive) BC_w -> f_w, BC_g -> f_g
	 * */
	
	/* boundary types */
	template<typename I> 
	std::vector< int >
	type( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
		  double time /*ndims*/,
		  double dt /*ndims*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		
		std::vector< int > bctype(Indices::numOfPVs, 0);
		
		
		if( property.mesh.isTopBoundary(globalPos)){
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet;
			bctype[indices.PVId_T] = indices.BCId_dirichlet;
			bctype[indices.PVId_C] = indices.BCId_dirichlet;
			bctype[indices.PVId_Sg] = indices.BCId_dirichlet;
		}
		// if( property.mesh.isBottomBoundary(globalPos)){
		// 	bctype[indices.PVId_Sg] = indices.BCId_neumann;
		// 	bctype[indices.PVId_Pw] = indices.BCId_neumann;
		// 	bctype[indices.PVId_T] = indices.BCId_neumann;
		// 	bctype[indices.PVId_C] = indices.BCId_neumann;
		// }
		// if( property.mesh.isLeftBoundary(globalPos)){
		// 	bctype[indices.PVId_Sg] = indices.BCId_neumann;
		// 	bctype[indices.PVId_Pw] = indices.BCId_neumann;
		// 	bctype[indices.PVId_T] = indices.BCId_neumann;
		// 	bctype[indices.PVId_C] = indices.BCId_neumann;
		// }
		// if( property.mesh.isWell(globalPos)){
		// 	bctype[indices.PVId_Pw] = indices.BCId_dirichlet;
		// }
		return bctype;
	}

	/* boundary values */
	template<typename I> 
	std::vector< double >
	value ( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
			double time/*s*/,
			double dt/*s*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		
		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();
		std::vector< double > bcvalue(Indices::numOfPVs,0.);
	    
		if( property.mesh.isTopBoundary(globalPos)){
			auto icv /*ndim*/ = icvalue.evaluate(cell_inside,iplocal);
			double Pw_top = icv[Indices::PVId_Pw] - (1000.* property.parameter.g()[dim-1] * press_rate * Xc_time * (time+dt)) / property.characteristicValue.P_c ;
							
			// std::cout << Pw_top << "  " << dt << std::endl;
			// exit(0);
			double Sg_top = icv[Indices::PVId_Sg];
			double xc_top = icv[Indices::PVId_C];
			double T_top  = icv[Indices::PVId_T] + 0.035 * press_rate * Xc_time * (time+dt) / property.characteristicValue.T_c;
			bcvalue[indices.PVId_Sg] = Sg_top;//property.parameter.InitialSg(globalPos);
			bcvalue[indices.PVId_Pw] = Pw_top ;//property.parameter.InitialPw(globalPos) + 1000 * 9.81 * 0.01 * Xc_time * (time+dt);//0.01 is the burial velocity m/year
			bcvalue[indices.PVId_T] = T_top;//property.parameter.InitialT(globalPos) ;//;//+ (2/(100*365*24*3600))*time - 3/2/1000 * globalPos[0]* time;
			bcvalue[indices.PVId_C] = xc_top;//property.parameter.InitialXC(globalPos);
		}
		
		if( property.mesh.isBottomBoundary(globalPos)){
			bcvalue[indices.PVId_T] = 0.035 * (property.characteristicValue.x_c/property.characteristicValue.T_c);
		}
		// if( property.mesh.isWell(globalPos)){
		// 	bcvalue[indices.PVId_Pw] = 8.e6/property.characteristicValue.P_c;
		// }

		return bcvalue;
	}

	// BC regarding phase Velocities
	template<typename I> 
	std::vector< int >
	velType( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
		  double time /*ndims*/,
		  double dt /*ndims*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		std::vector< int > bctype(Indices::numOfVelBCs, 0);

		if( property.mesh.isTopBoundary(globalPos)){
			bctype[indices.BCId_water] = indices.BCId_dirichlet;
			bctype[indices.BCId_heat] = indices.BCId_dirichlet;
			bctype[indices.BCId_salt] = indices.BCId_dirichlet;
			bctype[indices.BCId_gas] = indices.BCId_dirichlet;
		}
		// if( property.mesh.isBottomBoundary(globalPos)){
		// 	bctype[indices.BCId_water] = indices.BCId_neumann;
		// 	bctype[indices.BCId_heat] = indices.BCId_neumann;
		// 	bctype[indices.BCId_salt] = indices.BCId_neumann;
		// 	bctype[indices.BCId_gas] = indices.BCId_neumann;
		// }
		// if( property.mesh.isLeftBoundary(globalPos)){
		// 	bctype[indices.BCId_water] = indices.BCId_neumann;
		// 	bctype[indices.BCId_heat] = indices.BCId_neumann;
		// 	bctype[indices.BCId_salt] = indices.BCId_neumann;
		// 	bctype[indices.BCId_gas] = indices.BCId_neumann;
		// }
		// if( property.mesh.isWell(globalPos)){
		// 	bctype[indices.BCId_water] = indices.BCId_dirichlet;
		// }
		
		return bctype;
	}

	template<typename I> 
	std::vector< double >
	velValue ( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
			double time/*s*/,
			double dt/*s*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		
		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();

		std::vector< double > bcvalue(Indices::numOfVelBCs,0.);
		if( property.mesh.isTopBoundary(globalPos) ){
			auto icv /*ndim*/ = icvalue.evaluate(cell_inside,iplocal);

			double Pw_top = icv[Indices::PVId_Pw]//(property.parameter.InitialPw(globalPos)/property.characteristicValue.P_c ;
							- (1000. * property.parameter.g()[dim-1] * press_rate * Xc_time * (time+dt))/property.characteristicValue.P_c ;
			double Sg_top = icv[Indices::PVId_Sg];
			double xc_top = icv[Indices::PVId_C];
			double T_top  = icv[Indices::PVId_T] + 0.035 * press_rate * Xc_time * (time+dt) / property.characteristicValue.T_c;

			bcvalue[Indices::BCId_water] = Pw_top ;
			bcvalue[Indices::BCId_salt ] = xc_top ;
			bcvalue[Indices::BCId_heat ] = T_top  ;
			bcvalue[Indices::BCId_gas ] = Sg_top  ;

		}
		if( property.mesh.isBottomBoundary(globalPos) ){
			bcvalue[Indices::BCId_heat ] = 0.035 * (property.characteristicValue.x_c/property.characteristicValue.T_c);
		}
		// if( property.mesh.isWell(globalPos) ){
		// 	bcvalue[Indices::BCId_water ] = 8.e6/property.characteristicValue.P_c;
		// }
		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};


	

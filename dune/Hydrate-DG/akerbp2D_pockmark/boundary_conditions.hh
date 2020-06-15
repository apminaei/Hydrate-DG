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
	template<typename I> 
	std::vector< int >
	type( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
		  double time /*ndims*/,
		  double dt /*ndims*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		//std::cout << " iplocal_s = " << iplocal << " ip_global_s = " << globalPos << std::endl;
      //exit(0);
		std::vector< int > bctype(Indices::numOfPVs, 0);

		bctype[indices.PVId_Pc] = indices.BCId_dirichlet ;
		bctype[indices.PVId_Sh] = indices.BCId_dirichlet ;
		bctype[indices.PVId_T ] = indices.BCId_dirichlet ;
		bctype[indices.PVId_XCH4] = indices.BCId_dirichlet ;
		bctype[indices.PVId_YH2O ] = indices.BCId_dirichlet ;
		if( (time >= 0 ) and (time < (200*3600) )){
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		}
		else if ((time >= (200*3600) ) and (time < (350*3600)))
		{
			bctype[indices.PVId_Pw] = indices.BCId_neumann ;
		}
		else if (time >= (350*3600))	
		{
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		}
		// if(property.mesh.isLeftBoundary(globalPos)){
		// 	bctype[indices.PVId_Sg] = indices.BCId_dirichlet ;
		// 	bctype[indices.PVId_Sh] = indices.BCId_dirichlet ;
		// }
		
		// if( property.mesh.isBottomBoundary(globalPos) ){
			
		// 	bctype[indices.PVId_Sg] = indices.BCId_dirichlet ;
		// 	bctype[indices.PVId_Sh] = indices.BCId_dirichlet ;
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
		//const auto &cell_outside = intersection.outside();

		

		std::vector< double > bcvalue(Indices::numOfPVs,0.);
		
		bcvalue[indices.PVId_T ] = (4.+273.15)/ characteristicValues.T_c;//ProblemICValues(globalPos)[indices.PVId_T];//indices.BCId_neumann ;
		
		if( (time >= 0 ) and (time < (200*3600) )){
			bcvalue[indices.PVId_Pw] = 2.e6/characteristicValues.P_c ;
		}
		else if ((time >= (200*3600) ) and (time < (350*3600)))
		{
			//std::cout << " gas_gen = " << time << std::endl;
			
			bcvalue[indices.PVId_Pw] = 0. ;
		}
		else if ((time >= (350*3600) ) and (time < (450*3600)))	
		{
			bcvalue[indices.PVId_Pw] = 5.e6/characteristicValues.P_c ;
		}
		else if (time >= (450*3600)  )	
		{
			bcvalue[indices.PVId_Pw] = (5.e6 + 10.*(time- 450*3600))/characteristicValues.P_c ;
		}
		auto S = 0.0055;
		auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal);/*BrooksCParams[0] gives Pentry in Pa*/
		auto Sg = property.parameter.InitialSg(globalPos);
		auto Sh = property.parameter.InitialSh(globalPos);
		auto Sw = 1. - Sg - Sh;
		auto por = property.soil.SedimentPorosity(cell_inside, iplocal);
		double Pc = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal, Sw, Sh, por)
						* property.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4]);
		//auto Pg = property.parameter.InitialPw(globalPos) + Pc;
		//auto zCH4 = property.eos.EvaluateCompressibilityFactor(bcvalue[indices.PVId_T ] * characteristicValues.T_c, Pg * characteristicValues.P_c);
		//auto VLequil = property.mixture.EquilibriumMoleFractions( bcvalue[indices.PVId_T ]* characteristicValues.T_c, Pg * characteristicValues.P_c, S, zCH4);
		bcvalue[indices.PVId_XCH4] = property.parameter.InitialXCH4(globalPos);//VLequil[Indices::compId_XCH4];//indices.BCId_neumann ;
		bcvalue[indices.PVId_YH2O] = property.parameter.InitialYH2O(globalPos);//VLequil[Indices::compId_YH2O];
		bcvalue[indices.PVId_Pc] = 8.48e4 ;
		bcvalue[indices.PVId_Sh] = Sh ;
		// if(property.mesh.isLeftBoundary(globalPos)){
			
		// 	bcvalue[indices.PVId_Sg] = Sg ;
			
		// 	bcvalue[indices.PVId_Sh] = Sh ;
			
		// }
		
		// if( property.mesh.isBottomBoundary(globalPos) ){
			
		// 	bcvalue[indices.PVId_Sg] = Sg ;
		// 	bcvalue[indices.PVId_Sh] = Sh ;
		// }
		

		return bcvalue;
	}
	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }

};

// template<typename GV>
// class BCTypeParam0 : public Dune::PDELab::DirichletConstraintsParameters
// {
// private:
// 	const GV& gv ;
// 	//const Properties& property;
// 	const static int dim = GV::dimension;
// public:


//   //! construct from grid view
//   BCTypeParam0( const GV& gv_)
// 	: gv ( gv_ )
//   {
//   }

//   template<typename I>
//   bool isDirichlet(const I & intersection,
//                    const Dune::FieldVector<typename I::ctype, dim-1> & coord) const
//   {
//     Dune::FieldVector<typename I::ctype, dim>
//     x = intersection.geometry().global( coord );

//     //if( x[0] > 1. - 1.e-6 /*right*/ or x[1] < 0. + 1.e-6 /*bottom*/)
//         //return true; /*dirichlet*/
//     //else
//         //return false; /*left, top, bottom -> neumann*/
//  std::cout<< "boundary = " << x[0]  << ", " << x[1] << std::endl;
// 	return true;
//   }

// };
// template<typename GV>
// class BCTypeParam1 : public Dune::PDELab::DirichletConstraintsParameters
// {
// private:
// 	const GV& gv ;
// 	//const Properties& property;
// 	const static int dim = GV::dimension;
// public:


//   //! construct from grid view
//   BCTypeParam1(  const GV& gv_)
// 	: gv ( gv_ )
//   {
//   }

//   template<typename I>
//   bool isDirichlet(const I & intersection,
//                    const Dune::FieldVector<typename I::ctype, dim-1> & coord
//                    ) const
//   {
//     Dune::FieldVector<typename I::ctype, dim>
//     x = intersection.geometry().global( coord );

//     if( x[0] > 1. - 1.e-6 /*right*/ or x[1] < 0. + 1.e-6 )
//         return true; /*dirichlet*/
//     else
//         return false; /*left, top, bottom -> neumann*/

//   }

// };



	


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
	std::vector< int >
	type( const Dune::FieldVector<double,dim>& globalPos,
		  double time/*s*/,
		  double dt/*s*/ ) const {
		
		//auto xglobal = intersection.geometry().global(xlocal);
		
		std::vector< int > bctype(Indices::numOfPVs, 0);
		//std::vector< bool > bctypebool(Indices::numOfPVs, false);

		bctype[indices.PVId_T ] = indices.BCId_dirichlet ;
		//bctypebool[indices.PVId_T ] = true;
		//auto bcT = Dune::PDELab::makeBoundaryConditionFromCallable(gv);
		bctype[indices.PVId_C ] = indices.BCId_neumann ;
		//bctypebool[indices.PVId_C ] = false;
		bctype[indices.PVId_Sh] = indices.BCId_dirichlet ;
		//bctypebool[indices.PVId_Sh ] = false;
		bctype[indices.PVId_Pc] = indices.BCId_neumann ;
		//bctypebool[indices.PVId_Pc ] = false;
		bctype[indices.PVId_Sg] = indices.BCId_dirichlet ;
		//bctypebool[indices.PVId_Sg ] = false;
		bctype[indices.PVId_XCH4 ] = indices.BCId_dirichlet ;
		//bctypebool[indices.PVId_XCH4 ] = false;
		bctype[indices.PVId_YH2O ] = indices.BCId_dirichlet ;
		bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		//bctypebool[indices.PVId_YH2O ] = false;
		//if( property.mesh.isTopBoundary(globalPos) ){
		// if ((0. < time )  and (time <= (200.*3600))/* s */) {
		// 	bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		// 	//bctypebool[indices.PVId_Pw ] = true;
		// }
		// else if ( ((200.*3600) < time )  and (time <= (350.*3600)) /* s */) {
		// 	bctype[indices.PVId_Pw] = indices.BCId_neumann ;
		// 	//bctypebool[indices.PVId_Pw ] = false;
		// }
		// else {
		// 	bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
		// 	//bctypebool[indices.PVId_Pw ] = true;
		// }
		//  //std::cout<< "Pw_boundary = " << y << std::endl;
		
		return bctype;
	}

	/* boundary values */
	std::vector< double >
	value ( const Dune::FieldVector<double,dim>& globalPos,
			double time/*s*/,
			double dt/*s*/ ) const {

		//auto xglobal = intersection.geometry().global(xlocal);

		std::vector< double > bcvalue(Indices::numOfPVs,0.);
		
		bcvalue[indices.PVId_T ] = 4.+273.15 ;
		bcvalue[indices.PVId_Pw ] = 2.e6 ;
		bcvalue[indices.PVId_C ] = 0. ;
		bcvalue[indices.PVId_Sh] = 0.3 ;
		bcvalue[indices.PVId_Pc] = 8.48e4 ;
		bcvalue[indices.PVId_Sg] = 0. ;
		bcvalue[indices.PVId_XCH4 ] = 0. ;
		bcvalue[indices.PVId_YH2O ] = 0.0005 ;
		//if( property.mesh.isTopBoundary(globalPos) ){
		// if ((0. <= time )  and (time <= (200.*3600))/* second */) {
		// 	bcvalue[indices.PVId_Pw] = 2.e6 ; /* Pa */
		// }
		// else if ( ((200.*3600) < time )  and (time <= (350.*3600)) /* s */) {
		// 	bcvalue[indices.PVId_Pw] = 0. ;
		// }
		// else if ( ((350.*3600) < time )  and (time <= (450.*3600)) /* s */) {
		// 	bcvalue[indices.PVId_Pw] = 5.e6 ;
		// }
		// else{
		// 	bcvalue[indices.PVId_Pw] = 5.e6 + 10. * (time - 450. * 3600) ;
		// }
		
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
/*
* \brief constraints parameter class selecting boundary conditions
 */
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



	


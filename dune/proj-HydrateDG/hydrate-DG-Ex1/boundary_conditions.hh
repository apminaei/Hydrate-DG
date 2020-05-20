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
		
		if(property.mesh.isLeftBoundary(globalPos)){
			bctype[indices.PVId_T ] = indices.BCId_dirichlet ;
		
			bctype[indices.PVId_Sg] = indices.BCId_dirichlet ;
			
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
			bctype[indices.PVId_C] = indices.BCId_dirichlet ;
		}
		
		if( property.mesh.isRightBoundary(globalPos) ){
			bctype[indices.PVId_T ] = indices.BCId_dirichlet ;
		
			bctype[indices.PVId_Sg] = indices.BCId_dirichlet ;
			
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet ;
			bctype[indices.PVId_C] = indices.BCId_dirichlet ;
		}
		 //std::cout<< "Pw_boundary = " << y << std::endl;
		
		return bctype;
	}

	/* boundary values */
	std::vector< double >
	value ( const Dune::FieldVector<double,dim>& globalPos,
			double time/*s*/,
			double dt/*s*/ ) const {

		//auto xglobal = intersection.geometry().global(xlocal);

		std::vector< double > bcvalue(Indices::numOfPVs,0.);
		
		if(property.mesh.isLeftBoundary(globalPos)){
			bcvalue[indices.PVId_T ] = 60. + 273.15 ;
		
			bcvalue[indices.PVId_Sg] = 0.8 ;
			
			bcvalue[indices.PVId_Pw] = 1.5e6 ;
			bcvalue[indices.PVId_C] = 5.5e-3 ;
		}
		
		if( property.mesh.isRightBoundary(globalPos) ){
			bcvalue[indices.PVId_T ] = property.parameter.InitialT(globalPos) ;
		
			bcvalue[indices.PVId_Sg] = property.parameter.InitialSg(globalPos) ;
			
			bcvalue[indices.PVId_Pw] = property.parameter.InitialPw(globalPos); ;
			bcvalue[indices.PVId_C] = property.parameter.InitialXC(globalPos); ;
		}
		

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



	


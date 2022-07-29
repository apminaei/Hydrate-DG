/* ALL VALUES ARE NONDIMENSIONAL */
template <typename GV, typename Properties>
class ProblemBoundaryConditions
{
private:
	const GV &gv;
	const Properties &property;
	const static int dim = GV::dimension;
	Indices indices;
	CharacteristicValues characteristicValues;
	ProblemInitialConditions<GV, Properties> icvalue;
	double Xc_time = 1. / (36. * 24. * 36. * 1.e3);
	double press_rate = 1.e-3;

public:
	// ! construct from gridview
	ProblemBoundaryConditions(const GV &gv_, const Properties &property_)
		: gv(gv_),
		  property(property_),
		  icvalue(gv_, property_)
	{
	}
	/* boundary types */
	template <typename I>
	std::vector<int>
	type(I &intersection,
		 const Dune::FieldVector<double, dim - 1> &xlocal,
		 double time /*ndims*/,
		 double dt /*ndims*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		std::vector<int> bctype(Indices::numOfPVs, -1);

		if (property.mesh.isTopBoundary(globalPos))
		{
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet;
			bctype[indices.PVId_T] = indices.BCId_dirichlet;
			bctype[indices.PVId_C] = indices.BCId_dirichlet;
		}
		if (property.mesh.isBottomBoundary(globalPos))
		{

			bctype[indices.PVId_T] = indices.BCId_neumann;
			bctype[indices.PVId_C] = indices.BCId_neumann;
		}
		if (dim != 1)
		{
			if (property.mesh.isLeftBoundary(globalPos))
			{
				bctype[indices.PVId_T] = indices.BCId_neumann;
				bctype[indices.PVId_C] = indices.BCId_neumann;
			}
			if (property.mesh.isRightBoundary(globalPos))
			{
				bctype[indices.PVId_T] = indices.BCId_neumann;
				bctype[indices.PVId_C] = indices.BCId_neumann;
			}
		}

		return bctype;
	}

	/* boundary values */
	template <typename I>
	std::vector<double>
	value(I &intersection,
		  const Dune::FieldVector<double, dim - 1> &xlocal,
		  double time /*s*/,
		  double dt /*s*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();
		std::vector<double> bcvalue(Indices::numOfPVs, 0.);
		auto icv /*ndim*/ = icvalue.evaluate(cell_inside, iplocal);
		auto S = icv[Indices::PVId_C] * property.salt.MolarMass() / property.gas.MolarMass();
		auto T = icv[Indices::PVId_T] * property.characteristicValue.T_c;
		auto P = icv[Indices::PVId_Pw] * property.characteristicValue.P_c;
		if (property.mesh.isTopBoundary(globalPos))
		{

			double Pw_top = icv[Indices::PVId_Pw] + (920. * 9.81 * press_rate * Xc_time * (time + dt)) /*should increase */
														/ (property.characteristicValue.P_c);
			double T_top = icv[Indices::PVId_T] + (property.parameter.DTz() * press_rate * Xc_time * (time + dt)) / property.characteristicValue.T_c;

			double Sg_top = icv[Indices::PVId_Sg];
			double xc_top = icv[Indices::PVId_C];
			bcvalue[indices.PVId_Pw] = Pw_top;
			bcvalue[indices.PVId_T] = T_top;
			bcvalue[indices.PVId_C] = xc_top;
		}

		if (property.mesh.isBottomBoundary(globalPos))
		{
			bcvalue[indices.PVId_T] = -property.parameter.DTz() * (property.characteristicValue.x_c / property.characteristicValue.T_c); // with YASP should be negative
		}

		return bcvalue;
	}

	// BC regarding phase Velocities
	template <typename I>
	std::vector<int>
	velType(I &intersection, /*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double, dim - 1> &xlocal,
			double time /*ndims*/,
			double dt /*ndims*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		std::vector<int> bctype(Indices::numOfVelBCs, -1);

		if (property.mesh.isTopBoundary(globalPos))
		{
			// bctype[indices.BCId_water] = indices.BCId_dirichlet;
			bctype[indices.BCId_heat] = indices.BCId_dirichlet;
			bctype[indices.BCId_salt] = indices.BCId_dirichlet;
			// bctype[indices.BCId_gas] = indices.BCId_dirichlet;
		}
		if (property.mesh.isBottomBoundary(globalPos))
		{

			bctype[indices.BCId_water] = indices.BCId_neumann;
			bctype[indices.BCId_heat] = indices.BCId_neumann;
			bctype[indices.BCId_salt] = indices.BCId_neumann;
			bctype[indices.BCId_gas] = indices.BCId_neumann;
			// if( time > (50000. * 12.*30.*24*60.*60.)){
			// 	bctype[indices.BCId_heat] = indices.BCId_neumann;
			// }
		}
		if (dim != 1)
		{
			if (property.mesh.isLeftBoundary(globalPos))
			{
				bctype[indices.BCId_water] = indices.BCId_neumann;
				bctype[indices.BCId_heat] = indices.BCId_neumann;
				bctype[indices.BCId_salt] = indices.BCId_neumann;
				bctype[indices.BCId_gas] = indices.BCId_neumann;
			}
			if (property.mesh.isRightBoundary(globalPos))
			{
				bctype[indices.BCId_water] = indices.BCId_neumann;
				bctype[indices.BCId_heat] = indices.BCId_neumann;
				bctype[indices.BCId_salt] = indices.BCId_neumann;
				bctype[indices.BCId_gas] = indices.BCId_neumann;
			}
		}

		return bctype;
	}

	template <typename I>
	std::vector<double>
	velValue(I &intersection, /*const typename GV::Traits::template Codim<0>::Entity& element,*/
			 const Dune::FieldVector<double, dim - 1> &xlocal,
			 double time /*s*/,
			 double dt /*s*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();
		auto icv /*ndim*/ = icvalue.evaluate(cell_inside, iplocal);
		std::vector<double> bcvalue(Indices::numOfVelBCs, 0.);
		auto S = icv[Indices::PVId_C] * property.salt.MolarMass() / property.gas.MolarMass();
		auto T = icv[Indices::PVId_T] * property.characteristicValue.T_c;
		auto P = icv[Indices::PVId_Pw] * property.characteristicValue.P_c;
		if (property.mesh.isTopBoundary(globalPos))
		{

			// double Pw_top = icv[Indices::PVId_Pw];
			// double T_top  = icv[Indices::PVId_T];
			// if( time > (50000. * 12.*30.*24*60.*60.)){
			double Pw_top = icv[Indices::PVId_Pw] + (920. * 9.81 * press_rate * Xc_time * (time + dt)) /*should increase */
														/ (property.characteristicValue.P_c);
			double T_top = icv[Indices::PVId_T] + (property.parameter.DTz() * press_rate * Xc_time * (time + dt)) / property.characteristicValue.T_c;
			// }
			double Sg_top = icv[Indices::PVId_Sg]; // property.parameter.InitialSg(globalPos);
			double xc_top = icv[Indices::PVId_C];  // property.parameter.InitialXC(globalPos);

			bcvalue[Indices::BCId_water] = Pw_top;
			bcvalue[Indices::BCId_salt] = xc_top;
			bcvalue[Indices::BCId_heat] = T_top;
			// bcvalue[Indices::BCId_gas ] = Sg_top  ;
		}
		if (property.mesh.isBottomBoundary(globalPos))
		{

			// bcvalue[Indices::BCId_heat ] = icv[Indices::PVId_T];
			// if( time > (50000. * 12.*30.*24*60.*60.)){
			bcvalue[Indices::BCId_heat] = -property.parameter.DTz() * (property.characteristicValue.x_c / property.characteristicValue.T_c); // with YASP should be negative
																																			 // }
																																			 // bcvalue[Indices::BCId_water] = 0.;//icv[Indices::PVId_Pw];//- 0.06 * 1000. * property.parameter.g()[dim-1]/property.characteristicValue.P_c;
		}
		// if( property.mesh.isWell(globalPos) ){
		// 	bcvalue[Indices::BCId_water ] = 8.e6/property.characteristicValue.P_c;
		// }
		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV &getGridView() { return gv; }
};

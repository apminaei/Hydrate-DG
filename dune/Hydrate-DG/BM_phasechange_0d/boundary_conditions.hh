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
	double Xc_time = 100. * 3600.;

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
	type(I &intersection, const Dune::FieldVector<double, dim - 1> &xlocal, double time /* s */, double dt /* s */) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		std::vector<int> bctype(Indices::numOfPVs, -1);
		bctype[indices.PVId_T] = indices.BCId_dirichlet;
		bctype[indices.PVId_Sh] = indices.BCId_dirichlet;
		if ((time >= 0) and (time < (2.00 * Xc_time)))
		{
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet;
		}
		else if ((time >= (2.00 * Xc_time)) and (time < (3.50 * Xc_time)))
		{
			bctype[indices.PVId_Pw] = indices.BCId_neumann;
		}
		else if (time >= (3.50 * Xc_time))
		{
			bctype[indices.PVId_Pw] = indices.BCId_dirichlet;
		}

		return bctype;
	}

	/* boundary values */
	template <typename I>
	std::vector<double>
	value(I &intersection, const Dune::FieldVector<double, dim - 1> &xlocal, double time /*s*/, double dt /*s*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();
		std::vector<double> bcvalue(Indices::numOfPVs, 0.);

		bcvalue[indices.PVId_T] = property.parameter.InitialT(globalPos) / characteristicValues.T_c;

		if ((time >= 0) and (time < (2.00 * Xc_time)))
		{
			bcvalue[indices.PVId_Pw] = property.parameter.InitialPw(globalPos) / characteristicValues.P_c;
		}
		else if ((time >= (2.00 * Xc_time)) and (time < (3.50 * Xc_time)))
		{
			bcvalue[indices.PVId_Pw] = 0.;
		}
		else if ((time >= (3.50 * Xc_time)) and (time < (4.50 * Xc_time)))
		{
			bcvalue[indices.PVId_Pw] = 5.e6 / characteristicValues.P_c;
		}
		else if (time >= (4.50 * Xc_time))
		{
			bcvalue[indices.PVId_Pw] = (5.e6 + 10. * (time - 4.50 * Xc_time)) / characteristicValues.P_c;
		}

		auto Sh = property.parameter.InitialSh(globalPos);
		bcvalue[indices.PVId_Sh] = Sh;
		return bcvalue;
	}

	// BC regarding phase Velocities
	template <typename I>
	std::vector<int>
	velType(I &intersection, const Dune::FieldVector<double, dim - 1> &xlocal, double time /*s*/, double dt /*s*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		std::vector<int> bctype(Indices::numOfVelBCs, 0);
		bctype[indices.BCId_heat] = indices.BCId_dirichlet;
		if ((time >= 0) and (time < (2.00 * Xc_time)))
		{
			bctype[indices.BCId_water] = indices.BCId_dirichlet;
		}
		else if ((time >= (2.00 * Xc_time)) and (time < (3.50 * Xc_time)))
		{
			bctype[indices.BCId_water] = indices.BCId_neumann;
		}
		else if (time >= (3.50 * Xc_time))
		{
			bctype[indices.BCId_water] = indices.BCId_dirichlet;
		}

		return bctype;
	}

	template <typename I>
	std::vector<double>
	velValue(I &intersection, const Dune::FieldVector<double, dim - 1> &xlocal, double time /*s*/, double dt /*s*/) const
	{
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();

		std::vector<double> bcvalue(Indices::numOfVelBCs, 0.);
		bcvalue[indices.BCId_heat] = property.parameter.InitialT(globalPos) / characteristicValues.T_c;
		;
		if ((time >= 0) and (time < (2.00 * Xc_time)))
		{
			bcvalue[indices.BCId_water] = property.parameter.InitialPw(globalPos) / characteristicValues.P_c;
		}
		else if ((time >= (2.00 * Xc_time)) and (time < (3.50 * Xc_time)))
		{
			bcvalue[indices.BCId_water] = 0.;
		}
		else if ((time >= (3.50 * Xc_time)) and (time < (4.50 * Xc_time)))
		{
			bcvalue[indices.BCId_water] = 5.e6 / characteristicValues.P_c;
		}
		else if (time >= (4.50 * Xc_time))
		{
			bcvalue[indices.BCId_water] = (5.e6 + 10. * (time - 4.50 * Xc_time)) / characteristicValues.P_c;
		}
		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV &getGridView() { return gv; }
};

/*
 * Initial.hh
 * 
 *  class of initial primary variables that depends on 
 *  user defined initial values and properties of the problem
 *   
 *  All Values from initial_conditions are dimensional and must be transfered to nondim here 
 *    
 */

#ifndef INITIAL_HH_
#define INITIAL_HH_


/** \brief A function for initial values of Pw
 */
template<typename GV, typename Properties,  typename RF>
class Pw_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Pw_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Pw_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    y = icvalue.evaluate(e,xlocal)[Indices::PVId_Pw] ;	/* ndim */
    //std::cout<< "Pw_boundary = " << y << std::endl;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Pc
 */
// template<typename GV, typename Properties,  typename RF>
// class Pc_Initial
//   : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Pc_Initial<GV,Properties,RF> >
// {
// private:
// 	  const GV& gv;
//     const Properties& property;
//     ProblemInitialConditions<GV,Properties> icvalue;
// public:
//   typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

//   //! construct from grid view
//   Pc_Initial ( const GV& gv_, const Properties& property_ )
//   : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}

//   //! evaluate extended function on element
//   inline void evaluate (const typename Traits::ElementType& e,
//                         const typename Traits::DomainType& xlocal,
//                         typename Traits::RangeType& y) const
//   {

//     const int dim = Traits::GridViewType::Grid::dimension;
//     typedef typename Traits::GridViewType::Grid::ctype ctype;
//     Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

   
//     y = icvalue.evaluate(e,xlocal)[Indices::PVId_Pc]/CharacteristicValues::P_c ;	/* Pa */
//     //std::cout<< "Pc_boundary = " << y << std::endl;
//     return;
//   }
//   //! get a reference to the grid view
//   inline const GV& getGridView () {return gv;}
// };



/** \brief A function for initial values of Sg
 */
template<typename GV, typename Properties,  typename RF>
class Sg_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Sg_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Sg_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y= icvalue.evaluate(e,xlocal)[Indices::PVId_Sg] ; //initial water saturation
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};



/** \brief A function for initial values of Sh
 */
template<typename GV, typename Properties,  typename RF>
class Sh_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Sh_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Sh_Initial ( const GV& gv_ , const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_) {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y= icvalue.evaluate(e,xlocal)[Indices::PVId_Sh] ; //initial hydrate saturation

      		//std::cout << " Sh = " << y << std::endl;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of T
 */
template<typename GV, typename Properties,  typename RF>
class T_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, T_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  T_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y= icvalue.evaluate(e,xlocal)[Indices::PVId_T] ; //initial temperature
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/** \brief A function for initial values of XCH4
 */
template<typename GV, typename Properties,  typename RF>
class XCH4_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, XCH4_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  XCH4_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y= icvalue.evaluate(e,xlocal)[Indices::PVId_XCH4] ; //initial ch4 mole fraction /CharacteristicValues::x_c
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/** \brief A function for initial values of YH2O
 */
template<typename GV, typename Properties,  typename RF>
class YH2O_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, YH2O_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  YH2O_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y= icvalue.evaluate(e,xlocal)[Indices::PVId_YH2O] ; //initial h2o mole fraction /CharacteristicValues::x_c
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/** \brief A function for initial values of XC
 */
template<typename GV, typename Properties,  typename RF>
class XC_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, XC_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  XC_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y=icvalue.evaluate(e,xlocal)[Indices::PVId_C]; //initial salt mole fraction /CharacteristicValues::x_c 
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

#endif /* INITIAL_HH_ */

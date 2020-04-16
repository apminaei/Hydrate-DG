/*
 * Initial.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef INITIAL_HH_
#define INITIAL_HH_


/** \brief A function for initial values of Pg
 */
template<typename GV, typename RF>
class Pg_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Pg_Initial<GV,RF> >
{
private:
	  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Pg_Initial ( const GV& gv_ )
  : gv( gv_ ) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    IncludeClasses::ProblemSpecs problemSpecs;
    y = problemSpecs.ProblemICValues(x)[Indices::PVId_Pg]/CharacteristicValues::P_c ;	/* Pa */
    //std::cout<< "Pg_boundary = " << y << std::endl;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Pc
 */
template<typename GV, typename RF>
class Pc_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Pc_Initial<GV,RF> >
{
private:
	  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Pc_Initial ( const GV& gv_ )
  : gv( gv_ ) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    IncludeClasses::ProblemSpecs problemSpecs;
    y = problemSpecs.ProblemICValues(x)[Indices::PVId_Pc]/CharacteristicValues::P_c ;	/* Pa */
    //std::cout<< "Pc_boundary = " << y << std::endl;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};



/** \brief A function for initial values of Sw
 */
template<typename GV, typename RF>
class Sw_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Sw_Initial<GV,RF> >
{
private:
	  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Sw_Initial ( const GV& gv_ )
  : gv( gv_ ) {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    IncludeClasses::ProblemSpecs problemSpecs;
    y=problemSpecs.ProblemICValues(x)[Indices::PVId_Sw] ; //initial water saturation
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};



/** \brief A function for initial values of Sh
 */
template<typename GV, typename RF>
class Sh_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Sh_Initial<GV,RF> >
{
private:
	  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Sh_Initial ( const GV& gv_ )
  : gv( gv_ ) {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    IncludeClasses::ProblemSpecs problemSpecs;
    y=problemSpecs.ProblemICValues(x)[Indices::PVId_Sh] ; //initial hydrate saturation
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of T
 */
template<typename GV, typename RF>
class T_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, T_Initial<GV,RF> >
{
private:
	  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  T_Initial ( const GV& gv_ )
  : gv( gv_ ) {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    IncludeClasses::ProblemSpecs problemSpecs;
    y=problemSpecs.ProblemICValues(x)[Indices::PVId_T]/CharacteristicValues::T_c ; //initial temperature
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


#endif /* INITIAL_HH_ */

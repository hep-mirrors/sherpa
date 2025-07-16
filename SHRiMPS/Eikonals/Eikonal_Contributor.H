#ifndef SHRIMPS_Eikonals_Eikonal_Contributor_H
#define SHRIMPS_Eikonals_Eikonal_Contributor_H

#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Eikonals/Grid.H"
#include "SHRiMPS/Tools/DEQ_Solver.H"

namespace SHRIMPS {
  /*!
    \class Eikonal_Contributor
    \brief The representation of the two terms \f$\Omega_{i(k)}\f$ 
    contributing to the single channel eikonal.
  */
  class Eikonal_Contributor : public ATOOLS::Function_Base {
  private:
    Form_Factor      * p_ffi, * p_ffk;
    Eikonal_Parameters m_params;

    double m_b1, m_b2;
    double m_deltaff1, m_deltaff2, m_deltay;
    Grid   m_grid;

    bool Valid(const double & y) const;
  public:
    Eikonal_Contributor(Form_Factor * ff1,Form_Factor * ff2,
			const Eikonal_Parameters & params);
    ~Eikonal_Contributor() {}

    void   PrepareGrid(const int & ff1steps, const int & ff2steps);
    void   InsertValues(const size_t & i,const size_t & j,
			const std::vector<double> & values);
    void   SetB1B2(const double & b1,const double & b2);

    double operator()(const double & b1,const double & b2,const double & y);
    double operator()(const double & y);

    Form_Factor * FF1() const        { return p_ffi; } 
    Form_Factor * FF2() const        { return p_ffk; } 
    const double & Kappa() const     { return p_ffi->Kappa(); }
    const double & Lambda2() const   { return p_ffi->Lambda2(); }

    const double & B1max() const { return m_params.bmax; }
    const double & Y() const     { return m_params.originalY; }
  };
}

#endif

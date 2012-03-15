#ifndef SHRIMPS_Eikonals_Eikonal_Contributor_H
#define SHRIMPS_Eikonals_Eikonal_Contributor_H

#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Tools/DEQ_Solver.H"

namespace SHRIMPS {
  /*!
    \class Eikonal_Contributor
    \brief The representation of the two terms \f$\Omega_{i(k)}\f$ contributing to 
    the single channel eikonal.

    They are initialised and filled by the SHRIMPS::Eikonal_Creator with 
    help of the SHRIMPS::DEQ_Solver and then handed to the SHRIMPS::Omega_ik, who will take 
    care of deleting them.
  */
  class Eikonal_Contributor : public ATOOLS::Function_Base {
  private:
    Form_Factor * p_ff1, * p_ff2;
    double        m_originalY, m_cutoffY, m_Y;
    double        m_b1min, m_b2min, m_b1max, m_b2max, m_ff1max, m_ff2max;
    int           m_ff1steps, m_ff2steps, m_ysteps;
    double        m_deltaff1, m_deltaff2, m_deltay;
    double        m_b1, m_b2;
    std::vector<std::vector<std::vector<double> > > m_grid;

    double AnalyticalValue(const double & b1,const double & b2,const double & y);
  public:
    Eikonal_Contributor(Form_Factor * ff1,Form_Factor * ff2);
    ~Eikonal_Contributor();
    void   PrepareGrid(const int & ff1steps, const int & ff2steps);
    void   InsertValues(const int & i,const int & j,const std::vector<double> & values);
    double operator()(const double & b1,const double & b2,const double & y) {
      SetB1B2(b1,b2);
      return (*this)(y);
    }
    double operator()(const double & y) const;
    Form_Factor *  FF1() const   { return p_ff1; } 
    Form_Factor *  FF2() const   { return p_ff2; }
    const double & B1max() const { return m_b1max; }
    const double & Y() const     { return m_Y; }

    void  SetB1B2(const double & b1,const double & b2) { m_b1 = b1; m_b2 = b2; }
  };
}

#endif

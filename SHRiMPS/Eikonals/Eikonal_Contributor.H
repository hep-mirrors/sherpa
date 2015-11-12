#ifndef SHRIMPS_Eikonals_Eikonal_Contributor_H
#define SHRIMPS_Eikonals_Eikonal_Contributor_H

#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Tools/DEQ_Solver.H"

namespace SHRIMPS {
  /*!
    \class Eikonal_Contributor
    \brief The representation of the two terms \f$\Omega_{i(k)}\f$ 
    contributing to the single channel eikonal.
  */
  class Grid {
  private:
    size_t m_ff1steps, m_ff2steps;
    double m_ff1max, m_ff2max, m_Ymax;
    size_t m_ff1bin, m_ff2bin;
    double m_deltaff1, m_deltaff2, m_deltay;
    double m_ff1low, m_ff1up, m_ff2low, m_ff2up;
    double m_d1up,   m_d1low, m_d2up,   m_d2low;
    
    std::vector<std::vector<std::vector<double> > > m_grid;
    bool   FixBins(const double & ff1,const double & ff2);
    double ValueAtLowerYEdge();
    double ValueAtUpperYEdge();
    double Value(const double & y);
  public:
    Grid(const double & Ymax);
    ~Grid();
    
    void Initialise(const size_t & ff1steps,const size_t & ff2steps,
		    const double & ff1max,const double & ff2max);
    void InsertValues(const size_t & i,const size_t & j,
		      const std::vector<double> & values);
    double operator()(const double & ff1,const double & ff2,const double & y);
  };

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

    Form_Factor * FF() const         { return p_ffi; } 
    Form_Factor * FF1() const        { return p_ffi; } 
    const double & Kappa() const     { return p_ffi->Kappa(); }
    const double & Lambda2() const   { return p_ffi->Lambda2(); }



    Form_Factor *  FF2() const   { return 0; } //p_ff2; }
    const double & B1max() const { return m_params.bmax; }
    const double & Y() const     { return m_params.originalY; }
  };
}

#endif

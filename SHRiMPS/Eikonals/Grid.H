#ifndef SHRIMPS_Eikonals_Grid_H
#define SHRIMPS_Eikonals_Grid_H

#include "ATOOLS/Math/MathTools.H"
#include <vector>

namespace SHRIMPS {
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
}

#endif

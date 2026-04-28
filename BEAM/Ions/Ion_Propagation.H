#ifndef BEAM_Main_Ion_Propagation_H
#define BEAM_Main_Ion_Propagation_H

#include "BEAM/Ions/Ion_Potentials.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <vector>

namespace BEAM {
  struct Ion_Constituent;
  class  Ion_Base;
  
  class Ion_Propagation {
    Ion_Base     * p_ion;
    size_t         m_A, m_Z;
    double         m_R, m_V, m_rho0;
    double         m_sigma, m_sigma2, m_norm, m_rhonorm;
    Ion_Potentials m_potentials;

    size_t m_method, m_Nsteps;
    double m_deltat; 
    
    std::vector<Ion_Constituent *>  * p_constituents;
    std::vector<std::vector<ATOOLS::Vec3D>               > m_mom, m_pos;
    std::vector<std::vector<std::vector<ATOOLS::Vec3D> > > m_deltarvec;
    std::vector<std::vector<std::vector<double> > >        m_deltar,m_deltar2norm;
    std::vector<double> m_tilderho, m_rho, m_U;

    bool   m_analyse;
    std::map<std::string, ATOOLS::Histogram * > m_histos;
    void   InitAnalysis();
    void   FinishAnalysis();

    void   InitLists();
    void   InitParameters();
    void   UpdateDistances(size_t level=0);
    void   RungeKuttaStep(const double & dt);
    void   NewtonEulerStep(const double & dt);
    void   HarvestPositionsAndMomenta();
    void   CalculateRadii(int step=-1);
  public:
    Ion_Propagation(Ion_Base * ion);
    ~Ion_Propagation() { FinishAnalysis(); }
    
    void Initialise();
    void Evolve(int Nsteps=-1.,double d=-1.);
    void Step(double dt=-1.);

    const double  Pmax(size_t i) const;
    inline void   Update() {
      HarvestPositionsAndMomenta();
      UpdateDistances();
      m_potentials.UpdateDensities();
      m_potentials.UpdatePotentials();
    }
    inline const double & TildeRho(size_t i) const { return m_tilderho[i]; }
    inline const double &      Rho(size_t i) const { return m_tilderho[i]; }
  };
}

#endif
    

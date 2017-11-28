#ifndef ATOOLS_Phys_Particle_Dresser_H
#define ATOOLS_Phys_Particle_Dresser_H

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Selector_List.H"

namespace ATOOLS {
  class Particle_Dresser {
  private:
    bool m_on;
    Selector_List * p_sl;
    size_t m_algo;
    double m_exp, m_dR2global;
    std::map<kf_code,double> m_kfdR2s;
    std::vector<double> m_dR2;
    std::vector<size_t> m_photons;
    std::vector<size_t> m_charges;


    std::vector<double>               m_di,m_dj;
    std::vector<std::vector<double> > m_dij;

    void SetAlgo(std::string algo);

    double Pow(const double& x, const double& exp);
    double DeltaPhi(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    double DeltaR2(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);

    void ConeDress(Selector_List& sl);
    void RecombinationDress(Selector_List& sl);
  public:
    Particle_Dresser(const Flavour * fl,
                     const size_t& nin, const size_t& nout,
                     std::string algo, double dR);
    Particle_Dresser(std::string algo, double dR);
    ~Particle_Dresser();

    void SetFlavourDependentCone(kf_code kf,double dR) { m_kfdR2s[kf]=dR*dR; }

    void CompleteConeLists();

    ATOOLS::Vec4D_Vector Dress(const ATOOLS::Vec4D_Vector& p);
    void                 Dress(Selector_List& sl);

    inline bool IsOn() const { return m_on; }

    inline size_t Algorithm() const  { return m_algo; }
    inline double Exponent() const   { return m_exp; }
  };

}
#endif

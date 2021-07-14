#ifndef Recola_Virtual_H
#define Recola_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Process_Base.H"

namespace Recola {
  
  class Ioperator {
  private:
    bool m_on;
    size_t m_n;
    double m_mu2;
    std::vector<double> m_Q2i,m_mi,m_m2i,m_gi,m_Ki;
    std::vector<std::vector<double> > m_Gi;
    std::vector<std::vector<double> > m_Qik,m_VSik,m_VNSik;

    inline double Lambda(double x, double y, double z)
    { return x*x+y*y+z*z-2.*(x*y+x*z+y*z); }

    double VS(const size_t& i, const size_t& k,
              const double& sik, const double& q2ik,
              const size_t& id);
    double VNS(const size_t& i, const size_t& k,
              const double& sik, const double& q2ik);
    double Gamma(const size_t& i,const size_t& id);

  public:
    Ioperator(const ATOOLS::Flavour_Vector& fls, const bool &on);
    ~Ioperator();

    void ComputeConstants(const ATOOLS::Flavour_Vector& fls);
    double I(const ATOOLS::Vec4D_Vector& moms, const double& epsfac);

    inline void SetIRScale2(const double& mu2) { m_mu2=mu2; }
  };

  class Recola_Virtual : public PHASIC::Virtual_ME2_Base {

    unsigned int m_recola_id;
    bool m_ismapped;
    std::map<int,PHASIC::Process_Info> m_procmap;
    int m_modebackup;
    //double m_IRscale;
    //double m_UVscale;
    size_t m_voqcd, m_boqcd;
    Ioperator * p_iop;

  public:
 
    Recola_Virtual(const PHASIC::Process_Info& pi,
                   const ATOOLS::Flavour_Vector& flavs,
                   unsigned int recola_id);

    ~Recola_Virtual();

    inline bool SetColours(const ATOOLS::Vec4D_Vector& momenta) {
      return true;
    }


    inline double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

    bool IsMappableTo(const PHASIC::Process_Info& pi);
    
  };

}

#endif

#ifndef PHASIC_Channels_Channel_Elements_KK_H
#define PHASIC_Channels_Channel_Elements_KK_H

#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Selectors/Cut_Data.H"
namespace PHASIC {

  class Channel_Elements_KK {
    int    m_rs;
    double m_sran,m_weight;
  public:
    int     m_nin,m_nout,m_kkp,m_ed;
    double  m_r2,m_maxm2,m_maxn,m_mpss,m_gn,m_prevET,m_gam;
    /*!
      Determines the mass for the KK particle. This has to be done every time
      before a phase space point is generated.

      The distribution of the mass states is
      \f[
      m^{2}_{\vec{n}}=\frac{4\pi^{2}\vec{n}^{2}}{R^{2}},
      \f]
      where \f$\vec{n}\f$ is a $n$-dimensional vector of integers and with $R$ 
      contributes to the size of the extra dimensions, is given by
      \f[
      R^{n}=\frac{(4\pi)^{n/2}\Gamma(n/2)}{2M_s^{n+2}}.
      \f]
    */

    Channel_Elements_KK();
    ~Channel_Elements_KK();
    void Init(int nin,int nout,ATOOLS::Flavour* fl);
    void SetKKmass(double *ms, double ET, Cut_Data* cuts,double ran);
    double GetWeightKK(double& ran);

  };// end of class Channel_Elements_KK

  extern Channel_Elements_KK CEKK;                                             

}// end of namespace PHASIC     

#endif


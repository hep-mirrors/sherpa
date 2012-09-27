#ifndef PHASIC_Main_CS_Dipoles_H
#define PHASIC_Main_CS_Dipoles_H

#include "PHASIC++/Channels/CS_Dipole.H"

namespace PHASIC {

  class FF_Dipole: public CS_Dipole {
  private:
    double m_yexp, m_zexp;
  public:

    FF_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~FF_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

    void Calculate
    (const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
     const ATOOLS::Vec4D &pk,double &y,double &z,double &phi,
     ATOOLS::Vec4D &pijt,ATOOLS::Vec4D &pkt);
    void Construct
    (ATOOLS::Vec4D &pi,ATOOLS::Vec4D &pj,ATOOLS::Vec4D &pk,
     const double &y,const double &z,const double &phi,
     const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &pkt);

  };// end of class FF_Dipole

  class FI_Dipole: public CS_Dipole {
  private:
    double m_xmin, m_xexp, m_zexp;
  public:

    FI_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~FI_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

    void Calculate(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj,
                   const ATOOLS::Vec4D &pk,
                   double &x, double &z, double &phi,
                   ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt);
    void Construct(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
                   const double &x, const double &z, const double &phi,
                   const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt);

  };// end of class FI_Dipole

  class IF_Dipole: public CS_Dipole {
  private:
    double m_xmin,m_xexp, m_uexp;
    inline ATOOLS::Vec4D Rotate(const ATOOLS::Vec4D &p) const
    { return ATOOLS::Vec4D(p[0],-p[1],-p[2],-p[3]); }
  public:

    IF_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~IF_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

    void Calculate(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj,
                   const ATOOLS::Vec4D &pk,
                   double &x, double &u, double &phi,
                   ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt);
    void Construct(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
                   const double &x, const double &u, const double &phi,
                   const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt);

  };// end of class IF_Dipole

  class II_Dipole: public CS_Dipole {
  private:
    double m_xmin, m_xexp, m_vexp;
    inline ATOOLS::Vec4D Rotate(const ATOOLS::Vec4D &p) const
    { return m_ijt==m_sub.m_i?p:ATOOLS::Vec4D(p[0],-p[1],-p[2],-p[3]); }
  public:

    II_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~II_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

    void Calculate(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj,
                   const ATOOLS::Vec4D &pk,
                   const ATOOLS::Vec4D_Vector& kj,
                   double &x, double &v, double &phi,
                   ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt,
                   ATOOLS::Vec4D_Vector& kjt);
    void Construct(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
                   ATOOLS::Vec4D_Vector& kj,
                   const double &x, const double &v, const double &phi,
                   const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt,
                   const ATOOLS::Vec4D_Vector& kjt);

  };// end of class II_Dipole

}// end of namespace PHASIC

#endif

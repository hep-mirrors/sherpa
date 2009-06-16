#ifndef PHASIC_Main_CS_Dipoles_H
#define PHASIC_Main_CS_Dipoles_H

#include "PHASIC++/Channels/CS_Dipole.H"

namespace PHASIC {

  class FF_Dipole: public CS_Dipole {
  private:
    unsigned short int m_ymode, m_zmode;
    unsigned int m_exp;
  public:

    FF_Dipole(const ATOOLS::Cluster_Leg &lij,
              const ATOOLS::Cluster_Leg &lk,
              const ATOOLS::Cluster_Leg &li,
              const ATOOLS::Cluster_Leg &lj);
    ~FF_Dipole();

    CS_Dipole *Copy() const;

    bool GeneratePoint
    (ATOOLS::Cluster_Amplitude *const ampl,const double *rns);
    double GenerateWeight
    (const ATOOLS::Cluster_Amplitude *ampl);

    void Calculate
    (const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
     const ATOOLS::Cluster_Leg *lk,double &y,double &z,double &phi,
     ATOOLS::Vec4D &pijt,ATOOLS::Vec4D &pkt);
    void Construct
    (ATOOLS::Cluster_Leg *const li,ATOOLS::Cluster_Leg *const lj,
     ATOOLS::Cluster_Leg *const lk,ATOOLS::Cluster_Leg *const lij,
     const double &y,const double &z,const double &phi);

  };// end of class FF_Dipole


  class FI_Dipole: public CS_Dipole {
  private:
    unsigned short int m_xmode, m_zmode;
    unsigned int m_exp;
  public:

    FI_Dipole(const ATOOLS::Cluster_Leg &lij,
              const ATOOLS::Cluster_Leg &lk,
              const ATOOLS::Cluster_Leg &li,
              const ATOOLS::Cluster_Leg &lj);
    ~FI_Dipole();

    CS_Dipole *Copy() const;

    bool GeneratePoint
    (ATOOLS::Cluster_Amplitude *const ampl,const double *rns);
    double GenerateWeight
    (const ATOOLS::Cluster_Amplitude *ampl);

    void Calculate
    (const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
     const ATOOLS::Cluster_Leg *lk,double &y,double &z,double &phi,
     ATOOLS::Vec4D &pijt,ATOOLS::Vec4D &pkt);
    void Construct
    (ATOOLS::Cluster_Leg *const li,ATOOLS::Cluster_Leg *const lj,
     ATOOLS::Cluster_Leg *const lk,ATOOLS::Cluster_Leg *const lij,
     const double &y,const double &z,const double &phi);

  };// end of class FF_Dipole

  class IF_Dipole: public CS_Dipole {
  private:
    unsigned short int m_xmode, m_umode;
    unsigned int m_exp;
  public:

    IF_Dipole(const ATOOLS::Cluster_Leg &lij,
              const ATOOLS::Cluster_Leg &lk,
              const ATOOLS::Cluster_Leg &li,
              const ATOOLS::Cluster_Leg &lj);
    ~IF_Dipole();

    CS_Dipole *Copy() const;

    bool GeneratePoint
    (ATOOLS::Cluster_Amplitude *const ampl,const double *rns);
    double GenerateWeight
    (const ATOOLS::Cluster_Amplitude *ampl);

    void Calculate
    (const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
     const ATOOLS::Cluster_Leg *lk,double &y,double &z,double &phi,
     ATOOLS::Vec4D &pijt,ATOOLS::Vec4D &pkt);
    void Construct
    (ATOOLS::Cluster_Leg *const li,ATOOLS::Cluster_Leg *const lj,
     ATOOLS::Cluster_Leg *const lk,ATOOLS::Cluster_Leg *const lij,
     const double &y,const double &z,const double &phi);

  };// end of class FF_Dipole

  class II_Dipole: public CS_Dipole {
  private:
    unsigned short int m_xmode, m_vmode;
    unsigned int m_exp;
  public:

    II_Dipole(const ATOOLS::Cluster_Leg &lij,
              const ATOOLS::Cluster_Leg &lk,
              const ATOOLS::Cluster_Leg &li,
              const ATOOLS::Cluster_Leg &lj);
    ~II_Dipole();

    CS_Dipole *Copy() const;

    bool GeneratePoint
    (ATOOLS::Cluster_Amplitude *const ampl,const double *rns);
    double GenerateWeight
    (const ATOOLS::Cluster_Amplitude *ampl);

    void Calculate
    (const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
     const ATOOLS::Cluster_Leg *lk,double &y,double &z,double &phi,
     ATOOLS::Vec4D &pijt,ATOOLS::Vec4D &pkt);
    void Construct
    (ATOOLS::Cluster_Leg *const li,ATOOLS::Cluster_Leg *const lj,
     ATOOLS::Cluster_Leg *const lk,ATOOLS::Cluster_Leg *const lij,
     const double &y,const double &z,const double &phi);

  };// end of class FF_Dipole

}// end of namespace PHASIC

#endif

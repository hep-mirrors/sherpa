#ifndef SUSY_QCD_Splitting_Function_H
#define SUSY_QCD_Splitting_Function_H
//----------------------------------------------------------------------
//               QCD_Splitting_Functions.H
//----------------------------------------------------------------------
#include "Splitting_Function.H"

namespace APACIC {
  class Sudakov_Tools;

  // squark to squark + gluon  splitting function
  class SQuark__SQuark_Gluon : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
  public:
    SQuark__SQuark_Gluon(ATOOLS::Flavour squarkflavour);
    SQuark__SQuark_Gluon(ATOOLS::Flavour squarkflavour,Sudakov_Tools * tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double zmin, double zmax);
    double Integral(double zmin, double zmax);
  };

  // gluon to gluino + gluino  splitting function
  class Gluon__Gluino_Gluino : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
  public:
    Gluon__Gluino_Gluino();
    Gluon__Gluino_Gluino(Sudakov_Tools * tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double zmin, double zmax);
    double Integral(double zmin, double zmax);
  };

  // gluino to gluino + gluon  splitting function
  class Gluino__Gluino_Gluon : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
  public:
    Gluino__Gluino_Gluon();
    Gluino__Gluino_Gluon(Sudakov_Tools * tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double zmin, double zmax);
    double Integral(double zmin, double zmax);
  };

  // squark to squark + gluon  splitting function
  class Gluon__SQuark_SQuark : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
  public:
    Gluon__SQuark_SQuark(ATOOLS::Flavour squarkflavour);
    Gluon__SQuark_SQuark(ATOOLS::Flavour squarkflavour,Sudakov_Tools * tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double zmin, double zmax);
    double Integral(double zmin, double zmax);
  };

  /*!
    \file
    \brief contains the SYSY splitting fuctions 
  */


} // namespace APACIC
#endif









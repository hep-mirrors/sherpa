#ifndef APACIC_Showers_SUSY_QCD_Splitting_Functions_H
#define APACIC_Showers_SUSY_QCD_Splitting_Functions_H
//----------------------------------------------------------------------
//               QCD_Splitting_Functions.H
//----------------------------------------------------------------------
#include "APACIC++/Showers/Splitting_Function.H"

namespace APACIC {
  class Sudakov_Tools;

  // squark to squark + gluon  splitting function
  class SQuark__SQuark_Gluon : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
  public:
    SQuark__SQuark_Gluon(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour squarkflavour);
    SQuark__SQuark_Gluon(ATOOLS::Mass_Selector *&ms,
			 ATOOLS::Flavour squarkflavour,Sudakov_Tools * tools);
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
    Gluon__Gluino_Gluino(ATOOLS::Mass_Selector *&ms);
    Gluon__Gluino_Gluino(ATOOLS::Mass_Selector *&ms,Sudakov_Tools * tools);
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
    Gluino__Gluino_Gluon(ATOOLS::Mass_Selector *&ms);
    Gluino__Gluino_Gluon(ATOOLS::Mass_Selector *&ms,Sudakov_Tools * tools);
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
    Gluon__SQuark_SQuark(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour squarkflavour);
    Gluon__SQuark_SQuark(ATOOLS::Mass_Selector *&ms,
			 ATOOLS::Flavour squarkflavour,Sudakov_Tools * tools);
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









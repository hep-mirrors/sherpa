#ifndef Analysis_Observables_Four_Jet_Angles_H
#define Analysis_Observables_Four_Jet_Angles_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Four_Jet_Angle_Base : public Primitive_Observable_Base {
  public:

    Four_Jet_Angle_Base(unsigned int type,double xmin,double xmax,int nbins,
			const std::string & =std::string("Jets"));

    virtual double Calc(const std::vector<ATOOLS::Vec3D> & moms) = 0;
    void Evaluate(const ATOOLS::Blob_List & blobs,double value, double);
    void EvaluateNLOcontrib(double ,double);
    void EvaluateNLOevt();
  };

// 4JetBZ   93   0.    1.  20 Lin
// 4JetNR   93   0.    1.  20 Lin
// 4JetKSW  93   -1.   1.  20 Lin
// 4JetA34  93   -1.   1.  20 Lin

  class Bengtsson_Zerwas_Angle : public Four_Jet_Angle_Base {
  public:
    Bengtsson_Zerwas_Angle(unsigned int type,double xmin,double xmax,int nbins,
			   const std::string & =std::string("Jets"));

    double Calc(const std::vector<ATOOLS::Vec3D> & moms);
    Primitive_Observable_Base * Copy() const ;
  };

  class Nachtmann_Reiter_Angle : public Four_Jet_Angle_Base {
  public:
    Nachtmann_Reiter_Angle(unsigned int type,double xmin,double xmax,int nbins,
			   const std::string & =std::string("Jets"));

    double Calc(const std::vector<ATOOLS::Vec3D> & moms);
    Primitive_Observable_Base * Copy() const ;
  };

  class Koerner_Schierholz_Willrodt_Angle : public Four_Jet_Angle_Base {
  public:
    Koerner_Schierholz_Willrodt_Angle(unsigned int type,double xmin,double xmax,int nbins,
				      const std::string & =std::string("Jets"));

    double Calc(const std::vector<ATOOLS::Vec3D> & moms);
    Primitive_Observable_Base * Copy() const ;
  };

  class Alpha34_Angle : public Four_Jet_Angle_Base {
  public:
    Alpha34_Angle(unsigned int type,double xmin,double xmax,int nbins,
		  const std::string & =std::string("Jets"));

    double Calc(const std::vector<ATOOLS::Vec3D> & moms);
    Primitive_Observable_Base * Copy() const ;
  };

}

#endif

#ifndef Analysis_Observables_Jet_Cone_Distribution_H
#define Analysis_Observables_Jet_Cone_Distribution_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Triggers/Primitive_Calorimeter.H"
#include "AddOns/Analysis/Triggers/Calorimeter_Cone.H"
#include "ATOOLS/Math/Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Jet_Cone_Distribution : public Primitive_Observable_Base {
    double                           m_Etcut;
    double m_etamin, m_etamax;
    std::vector<Calorimeter_Cone *>  m_cones;
    std::vector<ATOOLS::Histogram *> m_histos;

    void   Fill(double weight=1.,double ncount=1);
  public:
    Jet_Cone_Distribution(const int linlog, const double Etcut, 
			  const double etamin, const double etamax,  
			  const double Rmin, const double Rmax, const int bins, 
			  Primitive_Analysis * const ana=0);
    ~Jet_Cone_Distribution();
    Primitive_Observable_Base * Copy() const ;
    void   Reset();
    void   Output(const std::string & pname);
    void   EndEvaluation(double scale=1.0);
    void   Restore(double scale=1.0);

    void   Evaluate(const ATOOLS::Particle_List &,double=1., double=1.0);
    void   Evaluate(const ATOOLS::Blob_List & blobs,double value=1., double=1.0);
  };

  inline void Jet_Cone_Distribution::Evaluate(const ATOOLS::Particle_List &,
					      double weight, double ncount) { Fill(weight,ncount); }
  inline void Jet_Cone_Distribution::Evaluate(const ATOOLS::Blob_List &,
					      double weight, double ncount) { Fill(weight,ncount); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class Jet_Cone_Dependence : public Primitive_Observable_Base {
    double                           m_Etcut;
    double                           m_etamin, m_etamax;
    int                              m_njetmin,m_njetmax;
    std::vector<Calorimeter_Cone *>  m_cones;
    std::vector<ATOOLS::Histogram *> m_histos;

    void   Fill(double weight=1.,double ncount=1);
  public:
    Jet_Cone_Dependence(const int linlog, const double Etcut, 
			const double etamin, const double etamax,  
			const double Rmin, const double Rmax, 
			const int njetmin, const int njetmax, 
			const int bins, 
			Primitive_Analysis * const);
    ~Jet_Cone_Dependence();
    Primitive_Observable_Base * Copy() const ;
    void   Reset();
    void   Output(const std::string & pname);
    void   EndEvaluation(double scale=1.0);
    void   Restore(double scale=1.0);

    void   Evaluate(const ATOOLS::Particle_List &,double=1., double=1.0);
    void   Evaluate(const ATOOLS::Blob_List & blobs,double value=1., double=1.0);
  };

  inline void Jet_Cone_Dependence::Evaluate(const ATOOLS::Particle_List &,
					    double weight, double ncount) { Fill(weight,ncount); }
  inline void Jet_Cone_Dependence::Evaluate(const ATOOLS::Blob_List &,
					    double weight, double ncount) { Fill(weight,ncount); }


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class Jet_Cone_Shape : public Primitive_Observable_Base {
    double                           m_Etcut, m_R;
    double                           m_etamin, m_etamax;
    int  m_jetmin, m_jetmax;
    Calorimeter_Cone               * p_cone;
    std::vector<ATOOLS::Histogram *> m_histos;
    void Fill(int jetno,double weight=1.,double ncount=1);
    void Fill(double weight=1.,double ncount=1);
  public:
    Jet_Cone_Shape(const int type,const double Etcut, const double radius, 
		   const double etamin, const double etamax,
		   const double Rmin, const double Rmax,
		   const int nbins, const int jetmin, const int jetmax, 
		   Primitive_Analysis *);
    ~Jet_Cone_Shape();
    Primitive_Observable_Base * Copy() const ;
    void   Reset();
    void   Output(const std::string & pname);
    void   EndEvaluation(double scale=1.0);
    void   Restore(double scale=1.0);

    void   Evaluate(const ATOOLS::Particle_List &,double=1., double=1.0);
    void   Evaluate(const ATOOLS::Blob_List & blobs,double value=1., double=1.0);
  };
  inline void Jet_Cone_Shape::Evaluate(const ATOOLS::Particle_List &,
				       double weight, double ncount) { Fill(weight,ncount); }
  inline void Jet_Cone_Shape::Evaluate(const ATOOLS::Blob_List &,
				       double weight, double ncount) { Fill(weight,ncount); }
}

#endif

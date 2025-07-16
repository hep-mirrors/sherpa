#ifndef Analysis_Observables_Jet_Mass_and_Broadening_H
#define Analysis_Observables_Jet_Mass_and_Broadening_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Observables/Momentum_Tensor.H"

namespace ANALYSIS {

  struct JetMass_Broadening_Data {
    friend std::ostream& operator<<( std::ostream&, const JetMass_Broadening_Data &);
    double heavyjetmass, lightjetmass;
    double widejetbroadening, narrowjetbroadening;
    JetMass_Broadening_Data(double mh, double ml,
		      double bw, double bn);
  };


  class JetMass_Broadening_Calculator : public Analysis_Object {
    std::string m_inkey, m_listname;
    std::string m_outkey;
  public:
    JetMass_Broadening_Calculator(const std::string & listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Analysis_Object * GetCopy() const;    
  };

  class Heavy_Jet_Mass : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Heavy_Jet_Mass(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Light_Jet_Mass : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Light_Jet_Mass(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Jet_Mass_Difference : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Jet_Mass_Difference(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Single_Inclusive_Jet_Mass : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Single_Inclusive_Jet_Mass(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Wide_Jet_Broadening : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Wide_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Narrow_Jet_Broadening : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Narrow_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Total_Jet_Broadening : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Total_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Jet_Broadening_Difference : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Jet_Broadening_Difference(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  class Single_Inclusive_Jet_Broadening : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Single_Inclusive_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

  inline JetMass_Broadening_Data::JetMass_Broadening_Data(double mh, double ml,
		      double bw, double bn)
    : heavyjetmass(mh), lightjetmass(ml), widejetbroadening(bw), narrowjetbroadening(bn)
  {
  }
  
} // namespace ANALYSIS
#endif

#ifndef Analysis_Observables_Sphericity_Observables_H
#define Analysis_Observables_Sphericity_Observables_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Observables/Momentum_Tensor.H"

namespace ANALYSIS {

  struct Sphericity_Data {
    friend std::ostream& operator<<( std::ostream&, const Sphericity_Data &);
    double        sphericity, aplanarity, planarity;
    Sphericity_Data(double spher,double aplan,double plan);
  };


  class Sphericity_Calculator : public Analysis_Object {
    Momentum_Tensor m_tensor;
    std::string m_key, m_listname;
  public:
    Sphericity_Calculator(const std::string & listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Analysis_Object * GetCopy() const;    
  };

  class Sphericity : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Sphericity(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class Aplanarity : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Aplanarity(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class Planarity : public Primitive_Observable_Base {
    std::string m_key;
  public:
    Planarity(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  inline  Sphericity_Data::Sphericity_Data(double spher,double aplan,double plan)
    : sphericity(spher), aplanarity(aplan), planarity(plan)
  {
  }
  
} // namespace ANALYSIS
#endif

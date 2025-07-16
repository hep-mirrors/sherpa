#ifndef Analysis_Observables_CD_Parameters_H
#define Analysis_Observables_CD_Parameters_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Observables/Momentum_Tensor.H"

namespace ANALYSIS {

  struct CD_Parameter_Data {
    friend std::ostream& operator<<( std::ostream&, const CD_Parameter_Data &);
    double        cparameter, dparameter;
    CD_Parameter_Data(double cpara,double dpara);
  };


  class CD_Parameter_Calculator : public Analysis_Object {
    Momentum_Tensor m_tensor;
    std::string m_key, m_listname;
  public:
    CD_Parameter_Calculator(const std::string & listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double value, double ncount);
    void EvaluateNLOevt();
    Analysis_Object * GetCopy() const;    
  };

  class C_Parameter : public Primitive_Observable_Base {
    std::string m_key;
  public:
    C_Parameter(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class D_Parameter : public Primitive_Observable_Base {
    std::string m_key;
  public:
    D_Parameter(int type, double xmin, double xmax, int nbin, std::string listname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };


  inline CD_Parameter_Data::CD_Parameter_Data(double cpara,double dpara)
    : cparameter(cpara), dparameter(dpara)
  {
  }
  
} // namespace ANALYSIS
#endif

#ifndef Analysis_Main_Analysis_Object_H
#define Analysis_Main_Analysis_Object_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Blob_List.H"
#include <string>

namespace ATOOLS {
  class Scoped_Settings;
}

namespace ANALYSIS {

  class Primitive_Analysis;

  struct Analysis_Key {
    Analysis_Key(ATOOLS::Scoped_Settings s, Primitive_Analysis* a):
      m_settings{s}, p_analysis{a}
    {}
    ATOOLS::Scoped_Settings m_settings;
    Primitive_Analysis *p_analysis;
  };

  class Analysis_Object {
  public:

    typedef ATOOLS::Getter_Function<Analysis_Object, Analysis_Key>
      Getter_Function;

  protected:

    Primitive_Analysis *p_ana;

    std::string m_name;
    bool        m_isobs, m_isdet;

  public:

    Analysis_Object();
    
    virtual ~Analysis_Object();

    virtual void Evaluate(const ATOOLS::Blob_List &blobs,
			  double value,double ncount) = 0;
    virtual void EvaluateNLOcontrib(double weight, double ncount);
    virtual void EvaluateNLOevt();

    virtual void Reset();
    virtual void Restore(double scale=1.0);
    virtual void EndEvaluation(double scale=1.0);
    virtual void Output(const std::string &pname);

    virtual void Test(const int mode=0);

    virtual Analysis_Object *GetCopy() const = 0;

    virtual Analysis_Object &operator+=(const Analysis_Object &obj);

    virtual void SetAnalysis(Primitive_Analysis *ana);
    virtual void SetName(const std::string &name);

    inline Primitive_Analysis *Analysis() const { return p_ana; }

    inline std::string Name() const { return m_name; }

    inline bool IsObservable() const { return m_isobs; }
    inline bool IsDetector()   const { return m_isdet; }

  };// end of class Analysis_Object

}// end of namespace ANALYSIS

#endif

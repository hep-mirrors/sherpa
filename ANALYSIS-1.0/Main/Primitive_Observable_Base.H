#ifndef Primitive_Observable_Base_H
#define Primitive_Observable_Base_H

#include "Selector.H"
#include "Vector.H"
#include "Flavour.H"
#include "Histogram.H"
#include "Getter_Function.H"

#include "Particle_List.H"
#include "Blob_List.H"

namespace ANALYSIS {

  class Primitive_Analysis;
  
  typedef std::vector<std::vector<std::string> > String_Matrix;

  class Primitive_Observable_Base {  
  public:

    typedef ATOOLS::Getter_Function<Primitive_Observable_Base,
				    String_Matrix> Getter_Function;

  protected:

    int                     m_type, m_nbins;
    double                  m_xmin, m_xmax;
    std::string             m_name, m_listname;
    ATOOLS::Histogram     * p_histo;
    int                     m_nout;
    ATOOLS::Flavour       * p_flavs;
    ATOOLS::Vec4D         * p_moms;
    std::string             m_blobtype;
    bool                    m_blobdisc, m_splitt_flag;

    Primitive_Analysis    * p_ana;
    ATOOLS::Selector_Base * p_sel;

  public:

    // constructors
    Primitive_Observable_Base();
    Primitive_Observable_Base(int _type,double _xmin,double _xmax,int _nbins,
			      ATOOLS::Selector_Base * _sel);
    Primitive_Observable_Base(Primitive_Observable_Base *);

    // destructor
    virtual ~Primitive_Observable_Base();

    // member functions
    virtual Primitive_Observable_Base *Copy() const = 0;

    virtual void Evaluate(int,const ATOOLS::Vec4D *,const ATOOLS::Flavour *,double=1.);
    virtual void Evaluate(const ATOOLS::Particle_List &,double=1.,int=1);
    virtual void Evaluate(const ATOOLS::Blob_List &,double=1.,int=1);

    virtual void EndEvaluation(double scale=1.);
    virtual void Reset();
    virtual void Output(const std::string & pname);

    virtual void SetBlobType(const std::string & =std::string(""));
    virtual void SetAnalysis(Primitive_Analysis*);

    virtual Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &);

    inline int    Type()  const { return m_type;  }
    inline int    Nbins() const { return m_nbins; }
    inline double Xmin()  const { return m_xmin;  }
    inline double Xmax()  const { return m_xmax;  }

    inline const std::string &Name() const       { return m_name;        }
    inline bool               Splittable() const { return m_splitt_flag; }

    inline ATOOLS::Histogram *const     Histo() const { return p_histo; }
    inline ATOOLS::Selector_Base *const Sel()   const { return p_sel;   }

    inline void SetName(const std::string &name) { m_name=name; }

  };// end of class Primitive_Observable_Base

}// end of namespace ANALYSIS

#endif

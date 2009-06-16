#ifndef Analysis_Observables_Primitive_Observable_Base_H
#define Analysis_Observables_Primitive_Observable_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Getter_Function.H"

#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Analysis/Main/Analysis_Object.H"

namespace ANALYSIS {

  int HistogramType(const std::string &scale);
  
  class Primitive_Observable_Base: public Analysis_Object {  
  public:

    typedef ATOOLS::Getter_Function<Primitive_Observable_Base,Argument_Matrix> 
    Observable_Getter_Function;

  protected:

    int                     m_type, m_nbins;
    double                  m_xmin, m_xmax;
    std::string             m_listname;
    ATOOLS::Histogram     * p_histo;
    int                     m_nout;
    ATOOLS::Flavour       * p_flavs;
    ATOOLS::Vec4D         * p_moms;
    std::string             m_blobtype;
    bool                    m_blobdisc, m_splitt_flag;

    mutable bool m_copied;

  public:

    // constructors
    Primitive_Observable_Base();
    Primitive_Observable_Base(int type,double xmin,double xmax,int nbins);
    Primitive_Observable_Base(const Primitive_Observable_Base & old);

    // destructor
    virtual ~Primitive_Observable_Base();

    // member functions
    virtual Primitive_Observable_Base *Copy() const = 0;

    Analysis_Object *GetCopy() const;

    virtual void Evaluate(int nout, const ATOOLS::Vec4D * moms,
			  const ATOOLS::Flavour * flavs,
			  double weight, double ncount);
    virtual void Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, double ncount);
    virtual void Evaluate(const ATOOLS::Blob_List & blobs, double value, double ncount);

    virtual void EndEvaluation(double scale=1.);
    virtual void Reset();
    virtual void Restore(double scale=1.0);
    virtual void Output(const std::string & pname);

    virtual void SetBlobType(const std::string & btype=std::string(""));
    virtual void SetAnalysis(Primitive_Analysis  * ana);

    virtual Primitive_Observable_Base & operator+=(const Primitive_Observable_Base & ob);
    virtual Analysis_Object & operator+=(const Analysis_Object & ob);

    int                     Type()  const;
    int                     Nbins() const;
    double                  Xmin()  const;
    double                  Xmax()  const;
    ATOOLS::Histogram     * Histo() const;
    bool                    Splittable() const;

    inline void SetName(const std::string &name) { m_name=name; }

  };// end of class Primitive_Observable_Base

  // --- inline functions ---
  inline int                     Primitive_Observable_Base::Type()       const { return m_type; }
  inline int                     Primitive_Observable_Base::Nbins()      const { return m_nbins; }
  inline double                  Primitive_Observable_Base::Xmin()       const { return m_xmin; }
  inline double                  Primitive_Observable_Base::Xmax()       const { return m_xmax; }

  inline ATOOLS::Histogram     * Primitive_Observable_Base::Histo()      const { return p_histo; }
  inline bool                    Primitive_Observable_Base::Splittable() const { return m_splitt_flag; }
}// end of namespace ANALYSIS

#endif

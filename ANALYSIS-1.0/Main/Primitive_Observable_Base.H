#ifndef Primitive_Observable_Base_H
#define Primitive_Observable_Base_H

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Selector.H"
#include "Vector.H"
#include "Flavour.H"
#include "Histogram.H"
#include "Message.H"

#include "Particle_List.H"
#include "Blob_List.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Primitive_Observable_Base {  
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
    Primitive_Observable_Base();
    Primitive_Observable_Base(int _type,double _xmin,double _xmax,int _nbins,
			      ATOOLS::Selector_Base * _sel);
    Primitive_Observable_Base(Primitive_Observable_Base *);
    virtual ~Primitive_Observable_Base();

    virtual Primitive_Observable_Base * Copy() const   =0; //   { return NULL; }

    virtual void Evaluate(int,const ATOOLS::Vec4D *,const ATOOLS::Flavour *,double=1.);
    virtual void Evaluate(const ATOOLS::Particle_List &,double=1.,int=1);
    virtual void Evaluate(const ATOOLS::Blob_List &,double=1.,int=1);

    virtual void EndEvaluation(double scale=1.);
    virtual void Reset();
    virtual void Output(const std::string & pname);

    //    virtual void SetFlavInfo(int,const Vec4D *,const Flavour *);
    virtual void SetBlobType(const std::string & =std::string(""));
    virtual void SetAnalysis(Primitive_Analysis*);

    virtual Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);

    int                     Type()  const;
    int                     Nbins() const;
    double                  Xmin()  const;
    double                  Xmax()  const;
    std::string             Name()  const;
    ATOOLS::Histogram     * Histo() const;
    ATOOLS::Selector_Base * Sel()   const;
    bool                    Splittable() const;
  };

  // --- inline functions ---
  inline int                     Primitive_Observable_Base::Type()       const { return m_type; }
  inline int                     Primitive_Observable_Base::Nbins()      const { return m_nbins; }
  inline double                  Primitive_Observable_Base::Xmin()       const { return m_xmin; }
  inline double                  Primitive_Observable_Base::Xmax()       const { return m_xmax; }
  inline std::string             Primitive_Observable_Base::Name()       const { return m_name; }
  inline ATOOLS::Histogram     * Primitive_Observable_Base::Histo()      const { return p_histo; }
  inline ATOOLS::Selector_Base * Primitive_Observable_Base::Sel()        const { return p_sel; }
  inline bool                    Primitive_Observable_Base::Splittable() const { return m_splitt_flag; }
}

#endif

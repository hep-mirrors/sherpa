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
  
  class String_Matrix: public std::vector<std::vector<std::string> > {
  private:

    Primitive_Analysis *p_analysis;

  public:

    // constructor
    inline String_Matrix() {}
    inline String_Matrix(const std::vector<std::vector<std::string> > &ref):
      std::vector<std::vector<std::string> >(ref) {}

    // member functions
    inline Primitive_Analysis *const operator()() const { return p_analysis; }
    inline const String_Matrix &operator()(Primitive_Analysis *const analysis) 
    { p_analysis=analysis; return *this; }
    
  };// end of class String_Matrix

  class Primitive_Observable_Base {  
  public:

    typedef ATOOLS::Getter_Function<Primitive_Observable_Base,String_Matrix> 
    Getter_Function;

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
    Primitive_Observable_Base(int type,double xmin,double xmax,int nbins,
			      ATOOLS::Selector_Base * sel=NULL);
    Primitive_Observable_Base(const Primitive_Observable_Base & old);

    // destructor
    virtual ~Primitive_Observable_Base();

    // member functions
    virtual Primitive_Observable_Base *Copy() const = 0;

    virtual void Evaluate(int nout, const ATOOLS::Vec4D * moms,
			  const ATOOLS::Flavour * flavs,
			  double weight, int ncount);
    virtual void Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
    virtual void Evaluate(const ATOOLS::Blob_List & blobs, double value, int ncount);

    virtual void EndEvaluation(double scale=1.);
    virtual void Reset();
    virtual void Output(const std::string & pname);

    virtual void SetBlobType(const std::string & btype=std::string(""));
    virtual void SetAnalysis(Primitive_Analysis  * ana);

    virtual Primitive_Observable_Base & operator+=(const Primitive_Observable_Base & ob);

    int                     Type()  const;
    int                     Nbins() const;
    double                  Xmin()  const;
    double                  Xmax()  const;
    std::string             Name()  const;
    void SetName(const std::string &);
    ATOOLS::Histogram     * Histo() const;
    ATOOLS::Selector_Base * Sel()   const;
    bool                    Splittable() const;

  };// end of class Primitive_Observable_Base

  // --- inline functions ---
  inline int                     Primitive_Observable_Base::Type()       const { return m_type; }
  inline int                     Primitive_Observable_Base::Nbins()      const { return m_nbins; }
  inline double                  Primitive_Observable_Base::Xmin()       const { return m_xmin; }
  inline double                  Primitive_Observable_Base::Xmax()       const { return m_xmax; }
  inline std::string             Primitive_Observable_Base::Name()       const { return m_name; }

  inline void Primitive_Observable_Base::SetName(const std::string & name ) { m_name=name; }

  inline ATOOLS::Histogram     * Primitive_Observable_Base::Histo()      const { return p_histo; }
  inline ATOOLS::Selector_Base * Primitive_Observable_Base::Sel()        const { return p_sel; }
  inline bool                    Primitive_Observable_Base::Splittable() const { return m_splitt_flag; }
}// end of namespace ANALYSIS

#endif

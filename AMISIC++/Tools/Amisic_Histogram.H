#ifndef AMISIC_Tools_Amisic_Histogram_H
#define AMISIC_Tools_Amisic_Histogram_H
/*!
  \file Amisic_Histogram.H
  \brief Declares the class Amisic_Histogram
*/

#include "ATOOLS/Math/Axis.H"

namespace AMISIC {

  struct hci {
  
    enum column {
      
      x_value  = 0,
      y_value  = 1,
      y_square = 2,
      maximum  = 3,
      entries  = 4,
      size     = 5
      
    };
    
  };// end of struct hci
  /*!
    \struct hci
    \brief Declares abbreviations for histogram column names.

    This struct is employed to declare abbreviations 
    for histogram column names.
  */

  template <class ArgumentType>
  class Amisic_Histogram {
  public:

    typedef ArgumentType Argument_Type;

    typedef std::vector<Argument_Type>   Argument_Vector;
    typedef std::vector<Argument_Vector> Argument_Matrix;
    
    typedef ATOOLS::Axis<Argument_Type> Axis_Type;

  private:

    size_t        m_nbins, m_extradim;
    Argument_Type m_xmin, m_xmax, m_entries;

    Argument_Matrix  m_data, m_sdata;
    Axis_Type       *p_xaxis, *p_yaxis;

    mutable Argument_Vector *p_integral;

    bool m_finished;

    std::string m_name;

    size_t FindX(const Argument_Type x) const;
    size_t FindY(const Argument_Type y) const;

    template <class Argument_Type> friend std::ostream &
    operator<<(std::ostream &str,
	       const Amisic_Histogram<Argument_Type> &histogram);

  public:

    // constructor
    Amisic_Histogram(const size_t extradim=0);

    // destructor
    ~Amisic_Histogram();

    // member functions
    bool Initialize(const Argument_Type xmin=0.0,
		    const Argument_Type xmax=0.0,
		    const size_t nbins=0);

    bool Import(const Argument_Matrix &ref,const bool overflow=false);
    bool Export(Argument_Matrix &ref,const bool overflow=false);

    size_t Add(Argument_Type value,
	       const Argument_Type weight,const size_t trials=1);
    size_t Set(Argument_Type value,
	       const Argument_Type weight,const size_t trials=1);

    Argument_Type Norm() const;

    void Scale(const Argument_Type scale);
    void ScaleExtra(const Argument_Type scale,const size_t dim=0);

    void Finish(); 
    void Clear();
      
    bool ReadIn(const std::string &filename,const std::string &m_datatag);
    bool WriteOut(const std::string &filename,const std::string &m_datatag,
		  const std::vector<std::string> &comments);

    Amisic_Histogram *GetIntegral(const bool reverse=false) const;
    Amisic_Histogram *GetDerivative() const;

    const Argument_Type operator()(const Argument_Type x) const;
    const Argument_Type operator[](const Argument_Type x) const;

    const Argument_Type GenerateX() const;

    void StoreData();
    void RestoreData();

    Argument_Type BinXMin(const size_t i); 
    Argument_Type BinXMax(const size_t i); 

    Argument_Type BinXMean(const size_t i); 
    Argument_Type BinVariance(const size_t i);
    Argument_Type BinError(const size_t i);

    Argument_Type BinContent(const size_t i);
    Argument_Type BinSumSqr(const size_t i);
    Argument_Type BinMax(const size_t i);
    Argument_Type BinEntries(const size_t i);
    Argument_Type BinExtra(const size_t i,const size_t dim=0);

    Argument_Type BinXMean(const Argument_Type &x); 
    Argument_Type BinVariance(const Argument_Type &x);
    Argument_Type BinError(const Argument_Type &x);

    Argument_Type BinContent(const Argument_Type &x);
    Argument_Type BinSumSqr(const Argument_Type &x);
    Argument_Type BinMax(const Argument_Type &x);
    Argument_Type BinEntries(const Argument_Type &x);
    Argument_Type BinExtra(const Argument_Type &x,const size_t dim=0);

    void SetBinContent(const size_t i,const Argument_Type &content);
    void SetBinSumSqr(const size_t i,const Argument_Type &sumsqr);
    void SetBinMax(const size_t i,const Argument_Type &max);
    void SetBinEntries(const size_t i,const Argument_Type &entries);
    void SetBinExtra(const size_t i,const Argument_Type &extra,
		     const size_t dim=0);

    void SetBinContent(const Argument_Type &x,const Argument_Type &content);
    void SetBinSumSqr(const Argument_Type &x,const Argument_Type &sumsqr);
    void SetBinMax(const Argument_Type &x,const Argument_Type &max);
    void SetBinEntries(const Argument_Type &x,const Argument_Type &entries);
    void SetBinExtra(const Argument_Type &x,const Argument_Type &extra,
		     const size_t dim=0);

    void AddBinContent(const size_t i,const Argument_Type &content);
    void AddBinSumSqr(const size_t i,const Argument_Type &sumsqr);
    void AddBinMax(const size_t i,const Argument_Type &max);
    void AddBinEntries(const size_t i,const Argument_Type &entries);
    void AddBinExtra(const size_t i,const Argument_Type &extra,
		     const size_t dim=0);

    void AddBinContent(const Argument_Type &x,const Argument_Type &content);
    void AddBinSumSqr(const Argument_Type &x,const Argument_Type &sumsqr);
    void AddBinMax(const Argument_Type &x,const Argument_Type &max);
    void AddBinEntries(const Argument_Type &x,const Argument_Type &entries);
    void AddBinExtra(const Argument_Type &x,const Argument_Type &extra,
		     const size_t dim=0);

    Axis_Type *XAxis() const;
    Axis_Type *YAxis() const;

    Argument_Type XMin() const;
    Argument_Type XMax() const;

    // inline functions
    size_t NBins() const;
    
    Argument_Type Entries() const;

    const std::string &Name() const;

    void SetXMin(const Argument_Type xmin);
    void SetXMax(const Argument_Type xmax);

    void SetNBins(const size_t nbins);

    void SetName(const std::string &name);

    void SetFinished(const bool finished);

  };// end of class Amisic_Histogram
  /*!
    \class Amisic_Histogram
    \brief Simple histogram class
    
    This class implements a histogram, which stores for each bin
    the sum of the weights \f$\sum w\f$, the sum of the squared
    weights \f$\sum w^2\f$, the maximum of the weights \f$max\{w\}\f$ 
    and the number of entries \f$n\f$.
    This allows for a determination of the mean, the variance and the 
    standard deviation for each bin.
  */

  template <class ArgumentType> std::ostream &
  operator<<(std::ostream &str,
	     const Amisic_Histogram<ArgumentType> &histogram);

}// end of namespace AMISIC

#endif

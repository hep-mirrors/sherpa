#ifndef Amisic_Histogram_H
#define Amisic_Histogram_H

#include "Axis.H"

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

  template <class ArgumentType>
  class Amisic_Histogram {
  public:

    typedef ArgumentType Argument_Type;

    typedef std::vector<Argument_Type>   Argument_Vector;
    typedef std::vector<Argument_Vector> Argument_Matrix;
    
    typedef ATOOLS::Axis<Argument_Type> Axis_Type;

  private:

    size_t        m_nbins;
    Argument_Type m_xmin, m_xmax, m_entries;

    Argument_Matrix  m_data;
    Axis_Type       *p_xaxis, *p_yaxis;

    bool m_finished;

    std::string m_name;

    size_t FindX(const Argument_Type x) const;
    size_t FindY(const Argument_Type y) const;

  public:

    // constructor
    Amisic_Histogram();

    // destructor
    ~Amisic_Histogram();

    // member functions
    bool Initialize(const Argument_Type xmin=0.0,
		    const Argument_Type xmax=0.0,
		    const size_t nbins=0);

    size_t Add(Argument_Type value,const Argument_Type weight);
    size_t Set(Argument_Type value,const Argument_Type weight);

    Argument_Type Norm() const;

    void Scale(const Argument_Type scale);

    void Finish(); 
    void Clear();
      
    bool ReadIn(const std::string &filename,const std::string &m_datatag);
    bool WriteOut(const std::string &filename,const std::string &m_datatag,
		  const std::vector<std::string> &comments);

    const Argument_Type operator()(const Argument_Type x) const;

    Argument_Type BinXMin(const size_t i); 
    Argument_Type BinXMax(const size_t i); 

    Argument_Type BinXMean(const size_t i); 
    Argument_Type BinError(const size_t i);

    Argument_Type BinContent(const size_t i);
    Argument_Type BinSumSqr(const size_t i);
    Argument_Type BinMax(const size_t i);
    Argument_Type BinEntries(const size_t i);

    Argument_Type BinContent(const Argument_Type x);
    Argument_Type BinSumSqr(const Argument_Type x);
    Argument_Type BinMax(const Argument_Type x);
    Argument_Type BinEntries(const Argument_Type x);

    void SetBinContent(const size_t i,const Argument_Type content);
    void SetBinSumSqr(const size_t i,const Argument_Type sumsqr);
    void SetBinMax(const size_t i,const Argument_Type max);
    void SetBinEntries(const size_t i,const Argument_Type entries);

    void SetBinContent(const Argument_Type x,const Argument_Type content);
    void SetBinSumSqr(const Argument_Type x,const Argument_Type sumsqr);
    void SetBinMax(const Argument_Type x,const Argument_Type max);
    void SetBinEntries(const Argument_Type x,const Argument_Type entries);

    // inline functions
    inline Argument_Type XMin() const { return m_xmin; }
    inline Argument_Type XMax() const { return m_xmax; }

    inline size_t NBins() const { return m_nbins; }

    inline Argument_Type Entries() const { return m_entries; }

    inline const std::string &Name() const { return m_name; }

    inline void SetXMin(const Argument_Type xmin) { m_xmin=xmin; }
    inline void SetXMax(const Argument_Type xmax) { m_xmax=xmax; }

    inline void SetNBins(const size_t nbins) { m_nbins=nbins; }

    inline void SetName(const std::string &name) { m_name=name; }

    inline void SetFinished(const bool finished) {m_finished=finished; }

    inline Axis_Type *const XAxis() const { return p_xaxis; }
    inline Axis_Type *const YAxis() const { return p_yaxis; }

  };// end of class Amisic_Histogram

}// end of namespace AMISIC

#endif

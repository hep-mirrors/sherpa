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

  public:

    // constructor
    Amisic_Histogram();

    // destructor
    ~Amisic_Histogram();

    // member functions
    void Initialize(const Argument_Type xmin=0.0,
		    const Argument_Type xmax=0.0,
		    const size_t nbins=0);

    size_t Add(Argument_Type value,const Argument_Type weight);
    size_t Set(Argument_Type value,const Argument_Type weight);

    Argument_Type Norm() const;

    void Finish(); 
      
    const Argument_Vector operator()(const Argument_Type x,
				     const hci::column col) const;

    // inline functions
    inline Argument_Type XMin() const { return m_xmin; }
    inline Argument_Type XMax() const { return m_xmax; }

    inline size_t NBins() const     { return m_nbins; }
    inline size_t Dimension() const { return m_dim;   }

    inline Axis_Type *const XAxis() const { return p_xaxis; }
    inline Axis_Type *const YAxis() const { return p_yaxis; }

    inline void SetXMin(const Argument_Type min) { m_xmin=xmin; }
    inline void SetXMax(const Argument_Type max) { m_xmax=xmax; }

    inline void SetNBins(const size_t nbins) { m_nbins=nbins; }

    inline Argument_Type BinXMin(const size_t i) 
    { return m_data[hci::x_value][i]; }
    inline Argument_Type BinXMax(const size_t i) 
    { return m_data[hci::x_value][i+1]; }

    inline Argument_Type BinContent(const size_t i)
    { return (*p_yaxis)[m_data[hci::y_value][i]]; }
    inline Argument_Type BinSumSqr(const size_t i)
    { return (*p_yaxis)[m_data[hci::y_square][i]]; }
    inline Argument_Type BinMax(const size_t i)
    { return (*p_yaxis)[m_data[hci::maximum][i]]; }
    inline Argument_Type BinEntries(const size_t i)
    { return (*p_yaxis)[m_data[hci::entries][i]]; }

    inline Argument_Type BinXMean(const size_t i) 
    { return (*p_xaxis)[((*p_xaxis)(m_data[hci::x_value][i+1])+
			 (*p_xaxis)(m_data[hci::x_value][i]))/2.0]; }
    inline Argument_Type BinError(const size_t i)
    { return (*p_yaxis)[m_data[hci::y_square][i]/
			m_data[hci::y_value][i]-
			m_data[hci::y_value][i]]; }

  };// end of class Amisic_Histogram

}// end of namespace AMISIC

#include "Amisic_Histogram.C"

#endif

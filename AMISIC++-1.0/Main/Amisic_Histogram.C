#ifndef Amisic_Histogram_C
#define Amisic_Histogram_C

#include "Amisic_Histogram.H"
#include "MathTools.H"
#include "My_Limits.H"

namespace AMISIC {
  
  template <class ArgumentType>
  Amisic_Histogram<ArgumentType>::Amisic_Histogram():
    m_entries(0.0),
    m_data(Argument_Matrix(hci::size)),
    p_xaxis(new Axis_Type()),
    p_yaxis(new Axis_Type()),
    m_finished(false) {}
  
  template <class ArgumentType>
  Amisic_Histogram<ArgumentType>::~Amisic_Histogram() 
  {
    delete p_yaxis;
    delete p_xaxis;
  }
  
  template <class ArgumentType>
  void Amisic_Histogram<ArgumentType>::Initialize(const Argument_Type xmin,
					   const Argument_Type xmax,
					   const size_t nbins)
  {
    if (nbins!=0) m_nbins=nbins;
    if (xmin!=xmax) {
      m_xmin=xmin;
      m_xmax=xmax;
    }
    for (size_t j=0;j<m_data.size();++j) m_data[j].resize(m_nbins+2);
    Argument_Type delta=((*p_xaxis)(m_xmax)-
			 (*p_xaxis)(m_xmin))/(double)nbins;
    for (size_t i=0;i<m_data[hci::x_value].size();++i) {
      for (size_t j=0;j<m_data.size();++j) m_data[j][i]=0.0;
      m_data[hci::x_value][i]=(*p_xaxis)[(*p_xaxis)(m_xmin)+(int(i)-1)*delta];
      m_data[hci::maximum][i]=-std::numeric_limits<Argument_Type>::max();
    }
    m_data[hci::x_value][0]=-std::numeric_limits<Argument_Type>::max();
    m_data[hci::x_value].back()=std::numeric_limits<Argument_Type>::max();
  }
  
  template <class ArgumentType>
  size_t Amisic_Histogram<ArgumentType>::Add(Argument_Type value,
					     const Argument_Type weight)
  {
    if (m_finished) return std::string::npos;
    ++m_entries;
    size_t l=0, r=m_data[hci::x_value].size()-1, i=(l+r)/2;
    double xi=m_data[hci::x_value][i];
    while (r-l>1) {
      if (value<xi) r=i;
      else l=i;
      i=(l+r)/2;
      xi=m_data[hci::x_value][i];
    }
    if (value<xi) --i;
    for (;i<m_data[hci::x_value].size();++i) {
      if (value>=m_data[hci::x_value][i] && 
	  value<m_data[hci::x_value][i+1]) {
	break;
      }
    }  
    m_data[hci::y_value][i]+=(*p_yaxis)(weight);
    m_data[hci::y_square][i]+=(*p_yaxis)(weight*weight);
    m_data[hci::maximum][i]=ATOOLS::Max(m_data[hci::maximum][i],
					(*p_yaxis)(weight));
    ++m_data[hci::entries][i];
    return i;
  }
  
  template <class ArgumentType>
  size_t Amisic_Histogram<ArgumentType>::Set(Argument_Type value,
					     const Argument_Type weight)
  {
    if (m_finished) return std::string::npos;
    ++m_entries;
    size_t l=0, r=m_data[hci::x_value].size()-1, i=(l+r)/2;
    double xi=m_data[hci::x_value][i];
    while (r-l>1) {
      if (value<xi) r=i;
      else l=i;
      i=(l+r)/2;
      xi=m_data[hci::x_value][i];
    }
    if (value<xi) --i;
    m_data[hci::y_value][i]=(*p_yaxis)(weight);
    m_data[hci::y_square][i]=(*p_yaxis)(weight*weight);
    m_data[hci::maximum][i]=(*p_yaxis)(weight);
    m_data[hci::entries][i]=1.;
    return i;
  }
  
  template <class ArgumentType>
  const std::vector<ArgumentType> 
  Amisic_Histogram<ArgumentType>::operator()(const Argument_Type x,
					     const hci::column col) const
  {
    x=(*p_xaxis)(x);
    size_t i=0;
    for (;i<m_data[0].size();++i) {
      if (x>=m_data[i] && x<m_data[i+1]) {
	break;
      }
    }  
    return (*p_yaxis)[m_data[col]];
  }
  
  template <class ArgumentType>
  ArgumentType Amisic_Histogram<ArgumentType>::Norm() const
  {
    Argument_Type integral=0.0;
    for (size_t i=0;i<m_data[hci::x_value].size()-1;++i) {
      integral+=((*p_xaxis)[m_data[hci::x_value][i+1]]-
		 (*p_xaxis)[m_data[hci::x_value][i]])*
	(*p_yaxis)[m_data[hci::yvalue][i]];
    }
    return integral;
  }
  
  template <class ArgumentType>
  void Amisic_Histogram<ArgumentType>::Finish()
  {
    for (size_t i=0;i<m_data[hci::x_value].size()-1;++i) {
      double binwidth=m_data[hci::x_value][i+1]-m_data[hci::x_value][i];
      m_data[hci::y_value][i]=(*p_yaxis)[m_data[hci::y_value][i]];
      m_data[hci::y_square][i]=(*p_yaxis)[m_data[hci::y_square][i]];
      m_data[hci::maximum][i]=(*p_yaxis)[m_data[hci::maximum][i]];
      m_data[hci::maximum][i]/=binwidth;
      m_data[hci::y_square][i]/=binwidth*m_entries;
      m_data[hci::y_value][i]/=binwidth*m_entries;
      m_data[hci::y_value][i]=(*p_yaxis)(m_data[hci::y_value][i]);
      m_data[hci::y_square][i]=(*p_yaxis)(m_data[hci::y_square][i]);
      m_data[hci::maximum][i]=(*p_yaxis)(m_data[hci::maximum][i]);
    }    
    m_finished=true;
  }

}// end of namespace AMISIC

#endif

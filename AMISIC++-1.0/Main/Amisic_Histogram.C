#include "Amisic_Histogram.H"
#include "Data_Writer.H"
#include "MathTools.H"
#include "My_Limits.H"

using namespace AMISIC;
  
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
void Amisic_Histogram<ArgumentType>::Clear() 
{
  m_entries=0.0;
  for (size_t i=0;i<hci::size;++i) {
    if (i!=hci::x_value) {
      m_data[i]=Argument_Vector(m_data[hci::x_value].size());
    }
  }
  m_finished=false;
}

template <class ArgumentType>
bool Amisic_Histogram<ArgumentType>::Initialize(const Argument_Type xmin,
						const Argument_Type xmax,
						const size_t nbins)
{
  if (nbins<1 || nbins>10000) return false;
  m_nbins=nbins;
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
  return true;
}

template <class ArgumentType>
size_t Amisic_Histogram<ArgumentType>::FindX(const Argument_Type x) const
{
  size_t l=0, r=m_data[hci::x_value].size()-1, i=(l+r)/2;
  double xi=m_data[hci::x_value][i];
  while (r-l>1) {
    if (x<xi) r=i;
    else l=i;
    i=(l+r)/2;
    xi=m_data[hci::x_value][i];
  }
  return l;
}

template <class ArgumentType>
size_t Amisic_Histogram<ArgumentType>::FindY(const Argument_Type y) const
{
  size_t l=0, r=m_data[hci::x_value].size()-1, i=(l+r)/2;
  double yi=m_data[hci::x_value][i];
  while (r-l>1) {
    if (y>yi) r=i;
    else l=i;
    i=(l+r)/2;
    yi=m_data[hci::x_value][i];
  }
  return l;
}

template <class ArgumentType>
size_t Amisic_Histogram<ArgumentType>::Add(Argument_Type value,
					   const Argument_Type weight)
{
  if (m_finished) return std::string::npos;
  ++m_entries;
  size_t i=FindX(value);
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
  size_t i=FindX(value);
  m_data[hci::y_value][i]=(*p_yaxis)(weight);
  m_data[hci::y_square][i]=(*p_yaxis)(weight*weight);
  m_data[hci::maximum][i]=(*p_yaxis)(weight);
  m_data[hci::entries][i]=1.;
  return i;
}
  
template <class ArgumentType>
const ArgumentType 
Amisic_Histogram<ArgumentType>::operator()(const Argument_Type x) const
{
  size_t l=FindX(x);
  if (l==0) ++l;
  else if (l+1==m_data[hci::x_value].size()-1) --l;
  double yl=m_data[hci::y_value][l];
  double ta=m_data[hci::y_value][l+1]-yl;
  double xl=(*p_xaxis)(m_data[hci::x_value][l]);
  ta/=(*p_xaxis)(m_data[hci::x_value][l+1])-xl;
  return (*p_yaxis)[yl+ta*((*p_yaxis)(x)-xl)];
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::Finish()
{
  if (m_finished) return;
  for (size_t i=0;i<m_data[hci::x_value].size();++i) {
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

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::Norm() const
{
  Argument_Type integral=0.0;
  for (size_t i=1;i<m_data[hci::x_value].size()-1;++i) {
    double width=m_data[hci::x_value][i+1]-m_data[hci::x_value][i];
    integral+=(*p_yaxis)[m_data[hci::y_value][i]]*width;
  }    
  return integral;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::Scale(const Argument_Type scale)
{
  for (size_t i=0;i<m_data[hci::y_value].size();++i) {
    m_data[hci::y_value][i]=(*p_yaxis)[m_data[hci::y_value][i]];
    m_data[hci::y_value][i]*=scale;
    m_data[hci::y_value][i]=(*p_yaxis)(m_data[hci::y_value][i]);
  }    
}

template <class ArgumentType>
bool Amisic_Histogram<ArgumentType>::ReadIn(const std::string &filename,
					    const std::string &datatag)
{
  if (filename=="") {
    ATOOLS::msg.Error()<<"Amisic_Histogram::ReadIn(..): "
		       <<"No filename specified. Abort."<<std::endl;
    return false;
  }
  std::string gridxscaling, gridyscaling, gridxvariable, gridyvariable;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","#");
  reader->SetInputFile(filename);
  reader->AddIgnore("!");
  reader->SetVectorType(reader->VHorizontal);
  if (!reader->ReadFromFile(m_xmin,"x_{min} :")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No x_{min} information in '"
		  <<filename<<"'."<<std::endl;
    return false;
  }
  if (!reader->ReadFromFile(m_xmax,"x_{max} :")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No x_{max} information in '"
		  <<filename<<"'."<<std::endl;
    return false;
  }
  if (!reader->ReadFromFile(gridxscaling,"x scale :")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No x scaling information in '"
		  <<filename<<"'."<<std::endl;
    return false;
  }
  if (!reader->ReadFromFile(gridyscaling,"y scale :")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No y scaling information in '"
		  <<filename<<"!."<<std::endl;
    return false;
  }
  std::vector<std::string> temp;
  if (!reader->VectorFromFile(temp,"x :")) gridxvariable="Unknown";
  else {
    gridxvariable=temp[0];
    for (unsigned int i=1;i<temp.size();++i) 
      gridxvariable+=std::string(" ")+temp[i];
  }
  if (!reader->VectorFromFile(temp,"y :")) gridyvariable="Unknown";
  else {
    gridyvariable=temp[0];
    for (unsigned int i=1;i<temp.size();++i) 
      gridyvariable+=std::string(" ")+temp[i];
  }
  p_xaxis->SetVariable(gridxvariable);
  p_yaxis->SetVariable(gridyvariable);
  p_xaxis->SetScaling(gridxscaling);
  p_yaxis->SetScaling(gridyscaling);
  reader->SetComment("!");
  reader->SetIgnore(":");
  std::vector<std::vector<double> > data;
  reader->MatrixFromFile(m_data,datatag);
  delete reader;
  if (m_data.size()<hci::size) {
    m_data.clear();
    return false;
  }
  m_nbins=m_data[hci::x_value].size()-2;
  m_entries=0.0;
  for (size_t i=0;i<m_nbins+2;++i) {
    m_entries+=m_data[hci::entries][i];
  }
  m_finished=true;
  return true;
}

template <class ArgumentType>
bool Amisic_Histogram<ArgumentType>::
WriteOut(const std::string &filename,const std::string &datatag,
	 const std::vector<std::string> &comments)
{
  if (filename=="") {
    msg_Tracking()<<"Amisic_Histogram::WriteOut(..): "
		  <<"No filename specified. Abort."<<std::endl;
    return false;
  }
  Finish();
  ATOOLS::Data_Writer *writer = new ATOOLS::Data_Writer(":",";","!");
  writer->SetOutputFile(filename);
  writer->SetBlank(ATOOLS::defaultblank);
  writer->WriteComment("===================="); 
  writer->WriteComment(" AMISIC++ grid file "); 
  writer->WriteComment("===================="); 
  writer->WriteComment(std::string("x scale : ")+
		       p_xaxis->Scaling()->Name());
  writer->WriteComment(std::string("y scale : ")+
		       p_yaxis->Scaling()->Name());
  writer->WriteComment("--------------------");
  writer->WriteComment(std::string("x : ")+
		       p_xaxis->Variable().Name());
  writer->WriteComment(std::string("y : ")+
		       p_yaxis->Variable().Name());
  writer->WriteComment("--------------------");
  writer->WriteComment(std::string("x_{min} : ")+
		       ATOOLS::ToString(m_xmin));
  writer->WriteComment(std::string("x_{max} : ")+
		       ATOOLS::ToString(m_xmax));
  writer->WriteComment("--------------------");
  writer->WriteComment(comments);
  writer->WriteComment("  Data Set follows  ");
  writer->WriteComment("--------------------");
  writer->SetBlank(ATOOLS::defaulttab);
  writer->MatrixToFile(m_data,datatag,true,
		       ATOOLS::nullstring,writer->MNormal,12);
  delete writer;
  return true;
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinXMin(const size_t i) 
{ 
  return m_data[hci::x_value][i]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinXMax(const size_t i) 
{ 
  return m_data[hci::x_value][i+1]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinContent(const size_t i)
{ 
  return (*p_yaxis)[m_data[hci::y_value][i]]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinSumSqr(const size_t i)
{ 
  return (*p_yaxis)[m_data[hci::y_square][i]]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinMax(const size_t i)
{ 
  return (*p_yaxis)[m_data[hci::maximum][i]]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinEntries(const size_t i)
{ 
  return (*p_yaxis)[m_data[hci::entries][i]]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinXMean(const size_t i) 
{ 
  return (*p_xaxis)[((*p_xaxis)(m_data[hci::x_value][i+1])+
		     (*p_xaxis)(m_data[hci::x_value][i]))/2.0]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinError(const size_t i)
{ 
  return (*p_yaxis)[m_data[hci::y_square][i]/
		    m_data[hci::y_value][i]-
		    m_data[hci::y_value][i]]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinContent(const Argument_Type x)
{ 
  return BinContent(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinSumSqr(const Argument_Type x)
{ 
  return BinSumSqr(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinMax(const Argument_Type x)
{ 
  return BinMax(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinEntries(const Argument_Type x)
{ 
  return BinEntries(FindX(x)); 
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinContent(const size_t i,const Argument_Type content)
{
  m_data[hci::y_value][i]=content;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinSumSqr(const size_t i,const Argument_Type sumsqr)
{
  m_data[hci::y_square][i]=sumsqr;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinMax(const size_t i,const Argument_Type max)
{
  m_data[hci::maximum][i]=max;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinEntries(const size_t i,const Argument_Type entries)
{
  m_data[hci::entries][i]=entries;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinContent(const Argument_Type x,const Argument_Type content)
{
  m_data[hci::y_value][FindX(x)]=content;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinSumSqr(const Argument_Type x,const Argument_Type sumsqr)
{
  m_data[hci::y_square][FindX(x)]=sumsqr;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinMax(const Argument_Type x,const Argument_Type max)
{
  m_data[hci::maximum][FindX(x)]=max;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinEntries(const Argument_Type x,const Argument_Type entries)
{
  m_data[hci::entries][FindX(x)]=entries;
}

template AMISIC::Amisic_Histogram<double>;

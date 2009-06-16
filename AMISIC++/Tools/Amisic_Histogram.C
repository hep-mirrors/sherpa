#include "AMISIC++/Tools/Amisic_Histogram.H"

#include "ATOOLS/Org/Data_Writer.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;
  
template <class ArgumentType> std::ostream &
AMISIC::operator<<(std::ostream &str,
		   const Amisic_Histogram<ArgumentType> &histogram) 
{
  str<<"("<<&histogram<<") {\n";
  {
    msg_Indent();
    str.precision(6);
    str<<std::setw(14)<<"value"
       <<std::setw(14)<<"weight"
       <<std::setw(14)<<"square weight"
       <<std::setw(14)<<"maximum"
       <<std::setw(14)<<"entries";
    for (size_t i=0;i<histogram.m_extradim;++i) 
      str<<std::setw(14)<<"extra "<<ATOOLS::ToString(i);
    str<<"\n";
    ATOOLS::Axis<ArgumentType> &yaxis=*histogram.YAxis();
    for (size_t i=0;i<histogram.m_data[hci::x_value].size();++i) {
      str<<std::setw(14)<<histogram.m_data[hci::x_value][i]<<" "
	 <<std::setw(14)<<yaxis[histogram.m_data[hci::y_value][i]]<<" "
	 <<std::setw(14)<<yaxis[histogram.m_data[hci::y_square][i]]<<" "
	 <<std::setw(14)<<yaxis[histogram.m_data[hci::maximum][i]]<<" "
	 <<std::setw(14)<<histogram.m_data[hci::entries][i];
      for (size_t j=0;j<histogram.m_extradim;++j) 
	str<<std::setw(14)<<" "<<histogram.m_data[hci::size+j][i];
      str<<"\n";
    }
  }
  return str<<"}"<<std::endl;
}

template <class ArgumentType>
Amisic_Histogram<ArgumentType>::Amisic_Histogram(const size_t extradim):
  m_extradim(extradim),
  m_entries(0.0),
  m_data(Argument_Matrix(hci::size+m_extradim)),
  p_xaxis(new Axis_Type()),
  p_yaxis(new Axis_Type()),
  p_integral(NULL),
  m_finished(false) {}

template <class ArgumentType>
Amisic_Histogram<ArgumentType>::~Amisic_Histogram() 
{
  if (p_integral!=NULL) delete p_integral;
  delete p_yaxis;
  delete p_xaxis;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::Clear() 
{
  m_entries=0.0;
  for (size_t i=0;i<m_data.size();++i) {
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
  if (nbins<2 || nbins>10000) return false;
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
  return true;
}

template <class ArgumentType>
bool Amisic_Histogram<ArgumentType>::Import(const Argument_Matrix &ref,
					     const bool overflow)
{
  if (ref.size()<2 || ref[0].size()<2 || ref[0].size() > 10000) return false;
  if (overflow) {
    m_nbins=ref[0].size()-2;
    m_xmin=ref[0][1];
    m_xmax=ref[0].back();
  }
  else {
    m_nbins=ref[0].size();
    m_xmin=ref[0][0];
    m_xmax=2.0*ref[0].back()-ref[0][m_nbins-2];
  }
  bool noverflow=!overflow;
  for (size_t j=0;j<m_data.size();++j) m_data[j].resize(m_nbins+2);
  for (size_t i=0;i<m_data[hci::x_value].size();++i) {
    for (size_t j=0;j<m_data.size();++j) m_data[j][i]=0.0;
    m_data[hci::maximum][i]=-std::numeric_limits<Argument_Type>::max();
    if (overflow^(i>0 && i<=ref[0].size())) {
      m_data[hci::x_value][i]=ref[0][i-noverflow];
      m_data[hci::y_value][i]=(*p_yaxis)(ref[1][i-noverflow]);
      if (ref.size()>2) 
	m_data[hci::y_square][i]=(*p_yaxis)(ref[2][i-noverflow]);
      if (ref.size()>3) 
	m_data[hci::maximum][i]=(*p_yaxis)(ref[3][i-noverflow]);
      if (ref.size()>4) m_data[hci::entries][i]=ref[4][i-noverflow];
    }
  }
  if (noverflow) {
    m_data[hci::x_value][0]=
      2.0*m_data[hci::x_value][1]-m_data[hci::x_value][2];
    m_data[hci::x_value].back()=
      2.0*m_data[hci::x_value][m_nbins]-m_data[hci::x_value][m_nbins-1];
  }
  return true;
}

template <class ArgumentType>
bool Amisic_Histogram<ArgumentType>::Export(Argument_Matrix &ref,
					     const bool overflow)
{
  if (ref.size()<2) return false;
  bool noverflow=!overflow;
  for (size_t j=0;j<ref.size();++j) ref[j].resize(m_nbins+2*overflow);
  for (size_t i=0;i<ref[0].size();++i) {
    ref[0][i]=m_data[hci::x_value][i+noverflow];
    ref[1][i]=(*p_yaxis)[m_data[hci::y_value][i+noverflow]];
    if (ref.size()>2) ref[2][i]=(*p_yaxis)[m_data[hci::y_square][i+noverflow]];
    if (ref.size()>3) ref[3][i]=(*p_yaxis)[m_data[hci::maximum][i+noverflow]];
    if (ref.size()>4) ref[4][i]=m_data[hci::entries][i+noverflow];
  }
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
					   const Argument_Type weight,
					   const size_t trials)
{
  if (m_finished) return std::string::npos;
  m_entries+=trials;
  size_t i=FindX(value);
  m_data[hci::y_value][i]+=(*p_yaxis)(weight);
  m_data[hci::y_square][i]+=(*p_yaxis)(weight*weight);
  m_data[hci::maximum][i]=ATOOLS::Max(m_data[hci::maximum][i],
				      (*p_yaxis)(weight));
  ++m_data[hci::entries][i];
  m_data[hci::entries][0]+=trials-1;
  return i;
}
  
template <class ArgumentType>
size_t Amisic_Histogram<ArgumentType>::Set(Argument_Type value,
					   const Argument_Type weight,
					   const size_t trials)
{
  if (m_finished) return std::string::npos;
  m_entries+=trials;
  size_t i=FindX(value);
  m_data[hci::y_value][i]=(*p_yaxis)(weight);
  m_data[hci::y_square][i]=(*p_yaxis)(weight*weight);
  m_data[hci::maximum][i]=(*p_yaxis)(weight);
  m_data[hci::entries][i]=1.;
  m_data[hci::entries][0]+=trials-1;
  return i;
}
  
template <class ArgumentType> const ArgumentType 
Amisic_Histogram<ArgumentType>::operator()(const Argument_Type x) const
{
  size_t l=FindX(x);
  if (l==0) ++l;
  else if (l+1==m_data[hci::x_value].size()-1) --l;
  double yl=m_data[hci::y_value][l];
  double ta=m_data[hci::y_value][l+1]-yl;
  double xl=(*p_xaxis)(m_data[hci::x_value][l]);
  ta/=(*p_xaxis)(m_data[hci::x_value][l+1])-xl;
  return (*p_yaxis)[yl+ta*((*p_xaxis)(x)-xl)];
}

template <class ArgumentType> const ArgumentType 
Amisic_Histogram<ArgumentType>::operator[](const Argument_Type y) const
{
  size_t l=0, r=m_data[hci::x_value].size()-1, i=(l+r)/2;
  bool inc=(*p_yaxis)[m_data[hci::y_value][l]]<
    (*p_yaxis)[m_data[hci::y_value][i]];
  double yi=m_data[hci::y_value][i], yc=(*p_yaxis)(y);
  while (r-l>1) {
    if (inc)
      if (yc<yi) r=i;
      else l=i;
    else 
      if (yc>yi) r=i;
      else l=i;
    i=(l+r)/2;
    yi=m_data[hci::y_value][i];
  }
  if (l==0) ++l;
  else if (l+1==m_data[hci::x_value].size()-1) --l;
  double yl=m_data[hci::y_value][l];
  double ta=m_data[hci::y_value][l+1]-yl;
  double xl=(*p_xaxis)(m_data[hci::x_value][l]);
  ta/=(*p_xaxis)(m_data[hci::x_value][l+1])-xl;
  return (*p_xaxis)[xl+(yc-yl)/ta];
}

template <class ArgumentType>
const ArgumentType Amisic_Histogram<ArgumentType>::DiceX() const
{
  if (p_integral==NULL) {
    p_integral = new Argument_Vector(m_nbins,0.0);
    double sum=0.0;
    for (size_t i=0;i<m_data[hci::x_value].size();++i) {
      double width=i<(m_data[hci::x_value].size()-1)?
	m_data[hci::x_value][i+1]-m_data[hci::x_value][i]:
	m_data[hci::x_value][i]-m_data[hci::x_value][i-1];
      (*p_integral)[i]=sum+=(*p_yaxis)[m_data[hci::y_value][i]]*width;
    }    
  }
  double y=ATOOLS::ran.Get()*(*p_integral)[m_data[hci::x_value].size()-1];
  size_t l=0, r=m_data[hci::x_value].size()-1, i=(l+r)/2;
  double yi=(*p_integral)[i];
  while (r-l>1) {
    if (y<yi) r=i;
    else l=i;
    i=(l+r)/2;
    yi=(*p_integral)[i];
  }
  if (l==0) ++l;
  else if (l+1==m_data[hci::x_value].size()-1) --l;
  double yl=(*p_integral)[l];
  double ta=(*p_integral)[l+1]-yl;
  double xl=(*p_xaxis)(m_data[hci::x_value][l]);
  ta/=(*p_xaxis)(m_data[hci::x_value][l+1])-xl;
  return (*p_xaxis)[xl+(y-yl)/ta];
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::Finish()
{
  if (m_finished) return;
  for (size_t i=0;i<m_data[hci::x_value].size();++i) {
    double binwidth=i<(m_data[hci::x_value].size()-1)?
      m_data[hci::x_value][i+1]-m_data[hci::x_value][i]:
      m_data[hci::x_value][i]-m_data[hci::x_value][i-1];
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
  for (size_t i=0;i<m_data[hci::x_value].size();++i) {
    double width=i<(m_data[hci::x_value].size()-1)?
      m_data[hci::x_value][i+1]-m_data[hci::x_value][i]:
      m_data[hci::x_value][i]-m_data[hci::x_value][i-1];
    integral+=(*p_yaxis)[m_data[hci::y_value][i]]*width;
  }    
  return integral;
}

template <class ArgumentType> Amisic_Histogram<ArgumentType> *const
Amisic_Histogram<ArgumentType>::GetIntegral(const bool reverse) const
{
  Amisic_Histogram<ArgumentType> *const integral = 
    new Amisic_Histogram<ArgumentType>();
  integral->p_xaxis->SetVariable(p_xaxis->Variable()->Name());
  integral->p_xaxis->SetScaling(p_xaxis->Scaling()->Name());
  integral->p_yaxis->SetVariable(p_yaxis->Variable()->Name());
  integral->p_yaxis->SetScaling(p_yaxis->Scaling()->Name());
  integral->Initialize(m_xmin,m_xmax,m_nbins);
  double sum=0.0;
  if (!reverse)
    for (size_t i=0;i<m_data[hci::x_value].size();++i) {
      double width=i<(m_data[hci::x_value].size()-1)?
	m_data[hci::x_value][i+1]-m_data[hci::x_value][i]:
	m_data[hci::x_value][i]-m_data[hci::x_value][i-1];
      double mean=m_data[hci::x_value][i]+width/2.0;
      integral->Set(mean,sum+=(*p_yaxis)[m_data[hci::y_value][i]]*width);
    }    
  else
    for (size_t i=m_data[hci::x_value].size();i>0;--i) {
      double width=i<m_data[hci::x_value].size()?
	m_data[hci::x_value][i]-m_data[hci::x_value][i-1]:
	m_data[hci::x_value][i-1]-m_data[hci::x_value][i-2];
      double mean=m_data[hci::x_value][i-1]+width/2.0;
      integral->Set(mean,sum+=(*p_yaxis)[m_data[hci::y_value][i-1]]*width);
    }    
  integral->SetFinished(true);
  return integral;
}

template <class ArgumentType> Amisic_Histogram<ArgumentType> *const
Amisic_Histogram<ArgumentType>::GetDerivative() const
{
  Amisic_Histogram<ArgumentType> *const derivative = 
    new Amisic_Histogram<ArgumentType>();
  derivative->p_xaxis->SetVariable(p_xaxis->Variable()->Name());
  derivative->p_xaxis->SetScaling(p_xaxis->Scaling()->Name());
  derivative->p_yaxis->SetVariable(p_yaxis->Variable()->Name());
  derivative->p_yaxis->SetScaling(p_yaxis->Scaling()->Name());
  derivative->Initialize(m_xmin,m_xmax,m_nbins);
  for (size_t i=1;i<m_data[hci::x_value].size()-1;++i) {
    double width=m_data[hci::x_value][i]-m_data[hci::x_value][i-1];
    double mean=m_data[hci::x_value][i]+width/2.0;
    derivative->Set(mean,((*p_yaxis)[m_data[hci::y_value][i]]-
			  (*p_yaxis)[m_data[hci::y_value][i-1]])/width);
  }    
  derivative->SetFinished(true);
  return derivative;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::Scale(const Argument_Type scale)
{
  for (size_t i=0;i<m_data[hci::y_value].size();++i) {
    m_data[hci::y_value][i]=(*p_yaxis)[m_data[hci::y_value][i]];
    m_data[hci::y_square][i]=(*p_yaxis)[m_data[hci::y_square][i]];
    m_data[hci::maximum][i]=(*p_yaxis)[m_data[hci::maximum][i]];
    m_data[hci::y_value][i]*=scale;
    m_data[hci::y_square][i]*=scale;
    m_data[hci::maximum][i]*=scale;
    m_data[hci::y_value][i]=(*p_yaxis)(m_data[hci::y_value][i]);
    m_data[hci::y_square][i]=(*p_yaxis)(m_data[hci::y_square][i]);
    m_data[hci::maximum][i]=(*p_yaxis)(m_data[hci::maximum][i]);
  }    
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
ScaleExtra(const Argument_Type scale,const size_t dim)
{
  if (dim>=m_extradim) return;
  for (size_t i=0;i<m_data[hci::y_value].size();++i) 
    m_data[hci::size+dim][i]*=scale;
}

template <class ArgumentType>
bool Amisic_Histogram<ArgumentType>::ReadIn(const std::string &filename,
					    const std::string &datatag)
{
  if (filename=="") {
    msg_Error()<<"Amisic_Histogram::ReadIn(..): "
		       <<"No filename specified. Abort."<<std::endl;
    return false;
  }
  std::string gridxscaling, gridyscaling, gridxvariable, gridyvariable;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader(" ",";",datatag,"=");
  reader->AddWordSeparator("\t");
  reader->SetAddCommandLine(false);
  reader->SetInputFile(filename);
  reader->AddIgnore("!");
  reader->AddIgnore(":");
  if (!reader->ReadFromFile(m_xmin,"x_{min}")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No x_{min} information in '"
		  <<filename<<"'."<<std::endl;
    return false;
  }
  if (!reader->ReadFromFile(m_xmax,"x_{max}")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No x_{max} information in '"
		  <<filename<<"'."<<std::endl;
    return false;
  }
  if (!reader->ReadFromFile(gridxscaling,"x-scale")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No x scaling information in '"
		  <<filename<<"'."<<std::endl;
    return false;
  }
  if (!reader->ReadFromFile(gridyscaling,"y-scale")) {
    msg_Tracking()<<"Amisic_Histogram::ReadIn(..): "
		  <<"No y scaling information in '"
		  <<filename<<"!."<<std::endl;
    return false;
  }
  std::vector<std::string> temp;
  if (!reader->VectorFromFile(temp,"x")) gridxvariable="Unknown";
  else {
    gridxvariable=temp[0];
    for (unsigned int i=1;i<temp.size();++i) 
      gridxvariable+=std::string(" ")+temp[i];
  }
  if (!reader->VectorFromFile(temp,"y")) gridyvariable="Unknown";
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
  reader->RescanInFile();
  std::vector<std::vector<double> > data;
  reader->SetMatrixType(ATOOLS::mtc::normal);
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
  ATOOLS::Data_Writer *writer = new ATOOLS::Data_Writer(" ",";","!",":");
  writer->SetOutputFile(filename);
  writer->SetBlank(ATOOLS::defaultblank);
  writer->WriteComment("===================="); 
  writer->WriteComment(" AMISIC++ grid file "); 
  writer->WriteComment("===================="); 
  writer->WriteComment(std::string("x-scale : ")+
		       p_xaxis->Scaling()->Name());
  writer->WriteComment(std::string("y-scale : ")+
		       p_yaxis->Scaling()->Name());
  writer->WriteComment("--------------------");
  writer->WriteComment(std::string("x : ")+
		       p_xaxis->Variable()->Name());
  writer->WriteComment(std::string("y : ")+
		       p_yaxis->Variable()->Name());
  writer->WriteComment("--------------------");
  writer->WriteComment(std::string("x_{min} : ")+
		       ATOOLS::ToString(m_xmin,12));
  writer->WriteComment(std::string("x_{max} : ")+
		       ATOOLS::ToString(m_xmax,12));
  writer->WriteComment("--------------------");
  writer->WriteComment(comments);
  writer->WriteComment("  Data Set follows  ");
  writer->WriteComment("--------------------");
  writer->SetBlank(ATOOLS::defaulttab);
  writer->MatrixToFile(m_data,datatag,true,
		       ATOOLS::nullstring,ATOOLS::mtc::normal,12);
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
  return m_data[hci::entries][i]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinXMean(const size_t i) 
{ 
  return (*p_xaxis)[((*p_xaxis)(m_data[hci::x_value][i+1])+
		     (*p_xaxis)(m_data[hci::x_value][i]))/2.0]; 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinVariance(const size_t i)
{ 
  return ((*p_yaxis)[m_data[hci::y_square][i]]-
	  ATOOLS::sqr((*p_yaxis)[m_data[hci::y_value][i]])/m_entries)/
    (m_entries-1.0); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinError(const size_t i)
{ 
  return sqrt(BinVariance(i)); 
}

template <class ArgumentType> ArgumentType Amisic_Histogram<ArgumentType>::
BinExtra(const size_t i,const size_t dim)
{
  return dim<m_extradim?m_data[hci::size+dim][i]:0.0;
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinXMean(const Argument_Type &x) 
{ 
  return BinXMean(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinVariance(const Argument_Type &x)
{ 
  return BinVariance(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinError(const Argument_Type &x)
{ 
  return BinError(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinContent(const Argument_Type &x)
{ 
  return BinContent(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinSumSqr(const Argument_Type &x)
{ 
  return BinSumSqr(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinMax(const Argument_Type &x)
{ 
  return BinMax(FindX(x)); 
}

template <class ArgumentType>
ArgumentType Amisic_Histogram<ArgumentType>::BinEntries(const Argument_Type &x)
{ 
  return BinEntries(FindX(x)); 
}

template <class ArgumentType> ArgumentType Amisic_Histogram<ArgumentType>::
BinExtra(const Argument_Type &x,const size_t dim)
{
  return BinExtra(FindX(x),dim); 
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinContent(const size_t i,const Argument_Type &content)
{
  m_data[hci::y_value][i]=(*p_yaxis)(content);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinSumSqr(const size_t i,const Argument_Type &sumsqr)
{
  m_data[hci::y_square][i]=(*p_yaxis)(sumsqr);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinMax(const size_t i,const Argument_Type &max)
{
  m_data[hci::maximum][i]=(*p_yaxis)(max);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinEntries(const size_t i,const Argument_Type &entries)
{
  m_data[hci::entries][i]=entries;
}

template <class ArgumentType> void Amisic_Histogram<ArgumentType>::
SetBinExtra(const size_t i,const Argument_Type &extra,const size_t dim)
{
  if (dim<m_extradim) m_data[hci::size+dim][i]=extra;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinContent(const Argument_Type &x,const Argument_Type &content)
{
  m_data[hci::y_value][FindX(x)]=(*p_yaxis)(content);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinSumSqr(const Argument_Type &x,const Argument_Type &sumsqr)
{
  m_data[hci::y_square][FindX(x)]=(*p_yaxis)(sumsqr);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinMax(const Argument_Type &x,const Argument_Type &max)
{
  m_data[hci::maximum][FindX(x)]=(*p_yaxis)(max);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
SetBinEntries(const Argument_Type &x,const Argument_Type &entries)
{
  m_data[hci::entries][FindX(x)]=entries;
}

template <class ArgumentType> void Amisic_Histogram<ArgumentType>::
SetBinExtra(const Argument_Type &x,const Argument_Type &extra,
	    const size_t dim)
{
  if (dim<m_extradim) m_data[hci::size+dim][FindX(x)]=extra;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinContent(const size_t i,const Argument_Type &content)
{
  m_data[hci::y_value][i]+=(*p_yaxis)(content);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinSumSqr(const size_t i,const Argument_Type &sumsqr)
{
  m_data[hci::y_square][i]+=(*p_yaxis)(sumsqr);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinMax(const size_t i,const Argument_Type &max)
{
  m_data[hci::maximum][i]+=(*p_yaxis)(max);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinEntries(const size_t i,const Argument_Type &entries)
{
  m_data[hci::entries][i]+=entries;
}

template <class ArgumentType> void Amisic_Histogram<ArgumentType>::
AddBinExtra(const size_t i,const Argument_Type &extra,const size_t dim)
{
  if (dim<m_extradim) m_data[hci::size+dim][i]+=extra;
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinContent(const Argument_Type &x,const Argument_Type &content)
{
  m_data[hci::y_value][FindX(x)]+=(*p_yaxis)(content);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinSumSqr(const Argument_Type &x,const Argument_Type &sumsqr)
{
  m_data[hci::y_square][FindX(x)]+=(*p_yaxis)(sumsqr);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinMax(const Argument_Type &x,const Argument_Type &max)
{
  m_data[hci::maximum][FindX(x)]+=(*p_yaxis)(max);
}

template <class ArgumentType>
void Amisic_Histogram<ArgumentType>::
AddBinEntries(const Argument_Type &x,const Argument_Type &entries)
{
  m_data[hci::entries][FindX(x)]+=entries;
}

template <class ArgumentType> void Amisic_Histogram<ArgumentType>::
AddBinExtra(const Argument_Type &x,const Argument_Type &extra,
	    const size_t dim)
{
  if (dim<m_extradim) m_data[hci::size+dim][FindX(x)]+=extra;
}

template <class ArgumentType> ATOOLS::Axis<ArgumentType> *const 
Amisic_Histogram<ArgumentType>::XAxis() const 
{ 
  return p_xaxis; 
}

template <class ArgumentType> ATOOLS::Axis<ArgumentType> *const 
Amisic_Histogram<ArgumentType>::YAxis() const 
{ 
  return p_yaxis; 
}

template <class ArgumentType> 
ArgumentType Amisic_Histogram<ArgumentType>::XMin() const 
{ 
  return m_xmin; 
}

template <class ArgumentType> 
ArgumentType Amisic_Histogram<ArgumentType>::XMax() const 
{ 
  return m_xmax; 
}

template <class ArgumentType> 
size_t Amisic_Histogram<ArgumentType>::NBins() const 
{ 
  return m_nbins; 
}
    
template <class ArgumentType> 
ArgumentType Amisic_Histogram<ArgumentType>::Entries() const 
{ 
  return m_entries; 
}

template <class ArgumentType> 
const std::string &Amisic_Histogram<ArgumentType>::Name() const 
{ 
  return m_name; 
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::SetXMin(const Argument_Type xmin) 
{ 
  m_xmin=xmin; 
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::SetXMax(const Argument_Type xmax) 
{
  m_xmax=xmax; 
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::SetNBins(const size_t nbins) 
{
  m_nbins=nbins; 
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::SetName(const std::string &name) 
{ 
  m_name=name; 
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::SetFinished(const bool finished) 
{ 
  m_finished=finished; 
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::StoreData()
{
  m_sdata=m_data;
}

template <class ArgumentType> 
void Amisic_Histogram<ArgumentType>::RestoreData()
{
  m_data=m_sdata;
  m_sdata=Argument_Matrix();
}

namespace AMISIC {

  template class AMISIC::Amisic_Histogram<double>;
  template std::ostream &
  operator<<<double>(std::ostream &str,
		     const Amisic_Histogram<double> &histogram);
  
}

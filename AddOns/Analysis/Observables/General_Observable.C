#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class General_Observable : public Primitive_Observable_Base {  
    std::string m_key;
  public:

    General_Observable(int type,double xmin,double xmax,int nbins,
		 const std::string & key);

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class General_Observable

}


#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(General_Observable_Getter,"Dummy",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base * 
General_Observable_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string ilist,olist,obskey;
  double keymin, keymax;
  int nbins;
  if (parameters.size()>0 && parameters[0].size()>=5) {
    obskey=parameters[0][0];
    keymin=ATOOLS::ToType<double>(parameters[0][1]);
    keymax=ATOOLS::ToType<double>(parameters[0][2]);
    nbins=ATOOLS::ToType<int>(parameters[0][3]);
    return new General_Observable(HistogramType(parameters[0][4]),
				  keymin,keymax,nbins,obskey);
  }
  return NULL;
}

void General_Observable_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"datakey  obsmin  obsmax nbins  Lin|LinErr|Log|LogErr"; 
}

General_Observable::General_Observable(int type,double xmin,double xmax,int nbins,
			   const std::string & key) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), m_key(key)
{
  
  m_name = m_key+".dat";
}
 
void General_Observable::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, double ncount)
{
  Blob_Data_Base * key=(*p_ana)[m_key];
  if (key) {
    double value = key->Get<double>();
    p_histo->Insert(value,weight,ncount); 
  }
  else {
    std::cout<<"warning #"<<m_key<<"# not found \n";
    p_histo->Insert(0.,0.,ncount); 
  }
}


Primitive_Observable_Base * General_Observable::Copy() const {
  return new General_Observable(m_type,m_xmin,m_xmax,m_nbins,m_key);
}



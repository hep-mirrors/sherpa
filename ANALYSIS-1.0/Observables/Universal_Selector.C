#include "Universal_Selector.H"
#include "MyStrStream.H"
#include "Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Universal_Selector_Getter,"UniSelector",
	       Primitive_Observable_Base,String_Matrix);

DECLARE_GETTER(General_Observable_Getter,"Dummy",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base * 
Universal_Selector_Getter::operator()(const String_Matrix &parameters) const
{
  std::string ilist,olist,obskey;
  double keymin, keymax;
  if (parameters.size()>0 && parameters[0].size()>=5) {
    obskey=parameters[0][0];
    keymin=ATOOLS::ToType<double>(parameters[0][1]);
    keymax=ATOOLS::ToType<double>(parameters[0][2]);
    ilist=parameters[0][3];
    olist=parameters[0][4];
    return new Universal_Selector(obskey,keymin,keymax,ilist,olist);
  }
  return NULL;
}

void Universal_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"obskey obsmin obsmax inlist outlist"; 
}



Universal_Selector::Universal_Selector(const std::string & obskey, const double keymin, const double keymax, 
				       const std::string & inlistname, const std::string & outlistname) :
  m_outlist(outlistname),m_key(obskey),m_keymin(keymin),m_keymax(keymax)
{
  m_name = std::string("UniversalSelector_")+outlistname;
  m_splitt_flag = false;
  m_listname = inlistname;
}

void Universal_Selector::CreateParticleList(bool force)
{
  Particle_List * pl_in = p_ana->GetParticleList(m_listname);
  if (pl_in==NULL) {
    msg.Out()<<"WARNING in Universal_Selector::Evaluate : particle list "<<m_listname<<" not found "<<std::endl;
  }
  Blob_Data_Base * key=(*p_ana)[m_key];

  Particle_List * pl_out = NULL;
  if (key && pl_in) {
    pl_out = new Particle_List;
    double value = key->Get<double>();
    if (m_keymin<=value && value<=m_keymax) {
      for (Particle_List::iterator pit=pl_in->begin();pit!=pl_in->end();++pit) {
	pl_out->push_back(new Particle(**pit));
      }
    }
  }
  if (!pl_out && force) pl_out = new Particle_List;
  if (pl_out) p_ana->AddParticleList(m_outlist,pl_out);
}

void Universal_Selector::Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout)
{
  CreateParticleList(false);
}

Primitive_Observable_Base *Universal_Selector::Copy() const
{
  return new Universal_Selector(m_key,m_keymin,m_keymax,m_listname,m_outlist);
}

// ======================================================================

Primitive_Observable_Base * 
General_Observable_Getter::operator()(const String_Matrix &parameters) const
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
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL), m_key(key)
{
  
  m_name = m_key+".dat";
}
 
void General_Observable::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, int ncount)
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



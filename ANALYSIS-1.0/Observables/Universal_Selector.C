#include "Universal_Selector.H"
#include "MyStrStream.H"
#include "Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Universal_Selector_Getter,"UniSelector",
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
  m_outlist(outlistname),m_key(obskey),m_keymin(keymin),m_keymax(keymin)
{
  std::cout<<"Universal_Selector::Universal_Selector("<<obskey<<","<<inlistname<<","<<outlistname<<")\n";

  m_splitt_flag = false;
  m_listname = inlistname;
}

void Universal_Selector::Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout)
{
  Particle_List * pl_in = p_ana->GetParticleList(m_listname);
  if (pl_in==NULL) {
    msg.Out()<<"WARNING in Universal_Selector::Evaluate : particle list "<<m_listname<<" not found "<<std::endl;
  }
  Particle_List * pl_out = new Particle_List;
  Blob_Data_Base * key=(*p_ana)[m_key];

  if (key && pl_in) {
    double value = key->Get<double>();
    if (m_keymin<=value && value<=m_keymax) {
      for (Particle_List::iterator pit=pl_in->begin();pit!=pl_in->end();++pit) {
	pl_out->push_back(new Particle(**pit));
      }
    }
  }
  p_ana->AddParticleList(m_outlist,pl_out);
}

Primitive_Observable_Base *Universal_Selector::Copy() const
{
  return new Universal_Selector(m_key,m_keymin,m_keymax,m_listname,m_outlist);
}

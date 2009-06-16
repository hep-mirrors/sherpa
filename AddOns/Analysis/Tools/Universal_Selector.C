#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "ATOOLS/Org/Message.H"

namespace ANALYSIS {

  class Universal_Selector: public Analysis_Object {  
  private:
    std::string m_listname, m_outlist;
    std::string m_key;
    double m_keymin, m_keymax;
  public:
    Universal_Selector(const std::string & obskey, const double keymin, const double keymax, 
		       const std::string & inlistname, const std::string & outlistname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight,double ncount);
    void CreateParticleList(bool force=true);
    Analysis_Object *GetCopy() const;

  };// end of class Universal_Selector

}// end of namespace ANALYSIS

#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Universal_Selector_Getter,"UniSelector",
	       Analysis_Object,Argument_Matrix);

Analysis_Object * 
Universal_Selector_Getter::operator()(const Argument_Matrix &parameters) const
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
  m_listname = inlistname;
}

void Universal_Selector::CreateParticleList(bool force)
{
  Particle_List * pl_in = p_ana->GetParticleList(m_listname);
  if (pl_in==NULL) {
    msg_Out()<<"WARNING in Universal_Selector::Evaluate : particle list "<<m_listname<<" not found "<<std::endl;
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

void Universal_Selector::Evaluate(const ATOOLS::Blob_List & ,double weight,double ncount)
{
  CreateParticleList(false);
}

Analysis_Object *Universal_Selector::GetCopy() const
{
  return new Universal_Selector(m_key,m_keymin,m_keymax,m_listname,m_outlist);
}


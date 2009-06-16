#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

namespace ANALYSIS {

  class Particle_Selector: public Analysis_Object {
  private:
    std::string  m_inlistname1, m_inlistname2, m_outlistname, m_qualname;
    int m_mode, m_init;
    ATOOLS::Particle_Qualifier_Base * p_qualifier;
  public:
    Particle_Selector(const std::string & inlistname1,
		      const std::string & inlistname2,
		      const std::string & outlistname,int mode);
    Particle_Selector(const std::string &inlist,const std::string &outlist,
		      const std::string &qualname);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List & ,double weight,double ncount);
    Analysis_Object *GetCopy() const;
    ~Particle_Selector();
  };

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(Particle_Selector_Getter,"PartSel",
 	       Analysis_Object,Argument_Matrix);

void Particle_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Qual    qualifier\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *
Particle_Selector_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string inlist="FinalState", outlist="Selected", qual("NotLepton");
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="Qual" && cur.size()>1) qual=cur[1];
  }
  return new Particle_Selector(inlist,outlist,qual);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

Particle_Selector::Particle_Selector(const std::string & inlistname1,
				     const std::string & inlistname2,
				     const std::string & outlistname, int mode)  : 
  m_inlistname1(inlistname1),m_inlistname2(inlistname2), 
  m_outlistname(outlistname), m_qualname(m_outlistname),
  m_mode(mode), m_init(0)
{
  m_name = std::string("ParticleSelector_")+outlistname;
  
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject(outlistname,"");
  if (!p_qualifier) {
    msg_Error()<<"ERROR in Particle_Selector: unknown particle qualifier '"
	       <<m_mode<<"'/'"<<outlistname<<"'"<<std::endl;
    p_qualifier = new Is_Charged();
  }
}

Particle_Selector::Particle_Selector(const std::string &inlist,
				     const std::string &outlist,
				     const std::string &qualname): 
  m_inlistname1(inlist), m_inlistname2(""), 
  m_outlistname(outlist), m_qualname(qualname),
  m_mode(0), m_init(1)
{
  m_name = outlist;
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::
    GetObject(m_qualname,m_qualname);
  if (p_qualifier==NULL) p_qualifier = new ATOOLS::Is_Not_Lepton(); 
}

void Particle_Selector::CreateParticleList()
{
  Particle_List * pl_in = NULL;
  if (m_mode<100) pl_in = p_ana->GetParticleList(m_inlistname1);
  else pl_in = p_ana->GetParticleList(m_inlistname2);
  if (pl_in==NULL) {
    msg_Out()<<"WARNING in Particle_Selector::Evaluate : particle list ";
    if (m_mode<100) msg_Out()<<m_inlistname1;
    else msg_Out()<<m_inlistname2;
    msg_Out()<<" not found "<<std::endl;
    return;
  }
  
  Particle_List * pl = new Particle_List;
  copy_if(pl_in->begin(),pl_in->end(),
	  back_inserter(*pl),*p_qualifier);
  
  p_ana->AddParticleList(m_outlistname,pl);
}

void Particle_Selector::Evaluate(const ATOOLS::Blob_List & ,double weight,double ncount)
{
  CreateParticleList();
}

Analysis_Object *Particle_Selector::GetCopy() const
{
  if (m_init==0) return new Particle_Selector
		   (m_inlistname1,m_inlistname2,m_outlistname,m_mode);
  return new Particle_Selector(m_inlistname1,m_outlistname,m_qualname);
}

Particle_Selector::~Particle_Selector()
{
  if (p_qualifier) delete p_qualifier;
}

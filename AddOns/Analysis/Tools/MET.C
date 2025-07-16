#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class MET_Define : public Analysis_Object {
  private:
    std::string  m_inlist, m_outlist;
    std::vector<Flavour> m_flavs;
  public:
    MET_Define(const std::string &inlist,
	       const std::string &outlist,
	       const std::vector<Flavour> &flavs);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List & ,double weight,double ncount);
    Analysis_Object *GetCopy() const;    
  };// end of class MET_Define

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(MET_Define,"METDefine",
	       Analysis_Object,Analysis_Key);

void ATOOLS::Getter<Analysis_Object,Analysis_Key,MET_Define>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList: list,\n"
     <<std::setw(width+7)<<" "<<"OutList: list,\n"
     <<std::setw(width+7)<<" "<<"Flavs: [flav1, .., flavN]\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *ATOOLS::Getter<Analysis_Object,Analysis_Key,MET_Define>::
operator()(const Analysis_Key& key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto inlist = s["InList"].SetDefault("FinalState").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("Selected").Get<std::string>();
  const auto rawflavs = s["Flavs"].SetDefault<int>({}).GetVector<int>();
  std::vector<Flavour> flavs;
  for (const auto& rawflav : rawflavs) {
    flavs.push_back(Flavour((kf_code)abs(rawflav)));
    if (rawflav<0) flavs.back()=flavs.back().Bar();
  }
  return new MET_Define(inlist,outlist,flavs);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

MET_Define::MET_Define(const std::string &inlist,
		       const std::string &outlist,
		       const std::vector<Flavour> &flavs):
  m_inlist(inlist), m_outlist(outlist), m_flavs(flavs) {}

void MET_Define::CreateParticleList()
{
  msg_Debugging()<<METHOD<<"() {\n";
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg_Error()<<METHOD<<"(): Missing lists: '"<<m_inlist
	       <<"'."<<std::endl;
    return;
  }
  Vec4D cms;
  Particle_List *outlist(new Particle_List());
  for (size_t i(0);i<inlist->size();++i) {
    bool found(false);
    for (size_t j(0);j<m_flavs.size();++j) {
      if (m_flavs[j].Includes((*inlist)[i]->Flav())) {
	msg_Debugging()<<"  skip "<<*(*inlist)[i]<<"\n";
	cms+=(*inlist)[i]->Momentum();
	found=true;
	break;
      }
    }
    if (!found) {
      msg_Debugging()<<"  keep "<<*(*inlist)[i]<<"\n";
      outlist->push_back(new Particle(*(*inlist)[i]));
    }
  }
  Particle *met(new Particle(1,Flavour(kf_none),cms));
  met->SetNumber(0);
  outlist->push_back(met);
  msg_Debugging()<<"  add  "<<*met<<"\n}\n";
  p_ana->AddParticleList(m_outlist,outlist);
}

void MET_Define::Evaluate(const ATOOLS::Blob_List & ,double weight,double ncount)
{
  CreateParticleList();
}

Analysis_Object *MET_Define::GetCopy() const
{
  return new MET_Define(m_inlist,m_outlist,m_flavs);
}


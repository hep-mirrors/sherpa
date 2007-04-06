#include "Analysis_Object.H"
#include "Particle_Qualifier.H"
#include "MyStrStream.H"
#include "Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Booster : public Analysis_Object {
  private:
    std::string  m_inlist, m_reflist, m_outlist;
    std::vector<Flavour> m_flavs;
    std::vector<int>     m_items;
  public:
    Booster(const std::string &inlist,
	    const std::string &reflist,
	    const std::string &outlist,
	    const std::vector<Flavour> &flavs,
	    const std::vector<int> &items);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout);
    Analysis_Object *GetCopy() const;    
  };// end of class Booster

  class Rotator : public Analysis_Object {
  private:
    std::string  m_inlist, m_reflist, m_outlist;
    std::vector<Flavour> m_flavs;
    std::vector<int>     m_items;
  public:
    Rotator(const std::string &inlist,
	    const std::string &reflist,
	    const std::string &outlist,
	    const std::vector<Flavour> &flavs,
	    const std::vector<int> &items);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout);
    Analysis_Object *GetCopy() const;    
  };// end of class Rotator

} // namespae ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(Booster_Getter,"CMSBoost",
 	       Analysis_Object,Argument_Matrix);

void Booster_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"RefList list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Flavs   flav1 .. flavN\n"
     <<std::setw(width+7)<<" "<<"Items   item1 .. itemN\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *
Booster_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Selected"), reflist("FinalState");
  std::vector<Flavour> flavs;
  std::vector<int> items;
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="RefList" && cur.size()>1) reflist=cur[1];
    else if (cur[0]=="Flavs" && cur.size()>1) {
      for (size_t i(1);i<cur.size();++i) {
	int kf(ToType<int>(cur[i]));
	flavs.push_back(Flavour((kf::code)abs(kf)));
	if (kf<0) flavs.back()=flavs.back().Bar();
      }
    }
    else if (cur[0]=="Items" && cur.size()>1) {
      for (size_t i(1);i<cur.size();++i) {
	items.push_back(ToType<int>(cur[i]));
      }
    }
  }
  items.resize(flavs.size(),0);
  return new Booster(inlist,reflist,outlist,flavs,items);
}

DECLARE_GETTER(Rotator_Getter,"ZRotate",
 	       Analysis_Object,Argument_Matrix);

void Rotator_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"RefList list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Flavs   flav1 .. flavN\n"
     <<std::setw(width+7)<<" "<<"Items   item1 .. itemN\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *
Rotator_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Selected"), reflist("FinalState");
  std::vector<Flavour> flavs;
  std::vector<int> items;
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="RefList" && cur.size()>1) reflist=cur[1];
    else if (cur[0]=="Flavs" && cur.size()>1) {
      for (size_t i(1);i<cur.size();++i) {
	int kf(ToType<int>(cur[i]));
	flavs.push_back(Flavour((kf::code)abs(kf)));
	if (kf<0) flavs.back()=flavs.back().Bar();
      }
    }
    else if (cur[0]=="Items" && cur.size()>1) {
      for (size_t i(1);i<cur.size();++i) {
	items.push_back(ToType<int>(cur[i]));
      }
    }
  }
  items.resize(flavs.size(),0);
  return new Rotator(inlist,reflist,outlist,flavs,items);
}

#include "Primitive_Analysis.H"

using namespace ATOOLS;

Booster::Booster(const std::string &inlist,
		 const std::string &reflist,
		 const std::string &outlist,
		 const std::vector<Flavour> &flavs,
		 const std::vector<int> &items):
  m_inlist(inlist), m_reflist(reflist), m_outlist(outlist), 
  m_flavs(flavs), m_items(items) {}

void Booster::CreateParticleList()
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  Particle_List *reflist(p_ana->GetParticleList(m_reflist));
  if (inlist==NULL || reflist==NULL) {
    msg.Error()<<METHOD<<"(): Missing lists: '"<<m_inlist
	       <<"','"<<m_reflist<<"'."<<std::endl;
    return;
  }
  Vec4D cms;
  for (size_t j(0);j<m_flavs.size();++j) {
    int no(-1);
    for (size_t i(0);i<reflist->size();++i) {
      if ((*reflist)[i]->Flav()==m_flavs[j]) {
	++no;
	if (no==m_items[j]) {
	  cms+=(*reflist)[i]->Momentum();
	  break;
	}
      }
    }
  }
  Poincare cmsboost(cms);
  Particle_List *outlist(new Particle_List());
  outlist->resize(inlist->size());
  for (size_t i(0);i<outlist->size();++i) {
    (*outlist)[i] = new Particle(*(*inlist)[i]);
    Vec4D p((*outlist)[i]->Momentum());
    cmsboost.Boost(p);
    (*outlist)[i]->SetMomentum(p);
  } 
  p_ana->AddParticleList(m_outlist,outlist);
}

void Booster::Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout)
{
  CreateParticleList();
}

Analysis_Object *Booster::GetCopy() const
{
  return new Booster(m_inlist,m_reflist,m_outlist,m_flavs,m_items);
}

Rotator::Rotator(const std::string &inlist,
		 const std::string &reflist,
		 const std::string &outlist,
		 const std::vector<Flavour> &flavs,
		 const std::vector<int> &items):
  m_inlist(inlist), m_reflist(reflist), m_outlist(outlist), 
  m_flavs(flavs), m_items(items) {}

void Rotator::CreateParticleList()
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  Particle_List *reflist(p_ana->GetParticleList(m_reflist));
  if (inlist==NULL || reflist==NULL) {
    msg.Error()<<METHOD<<"(): Missing lists: '"<<m_inlist
	       <<"','"<<m_reflist<<"'."<<std::endl;
    return;
  }
  Vec4D zaxis;
  for (size_t j(0);j<m_flavs.size();++j) {
    int no(-1);
    for (size_t i(0);i<reflist->size();++i) {
      if ((*reflist)[i]->Flav()==m_flavs[j]) {
	++no;
	if (no==m_items[j]) {
	  zaxis+=(*reflist)[i]->Momentum();
	  break;
	}
      }
    }
  }
  Poincare zrot(zaxis,Vec4D::ZVEC);
  Particle_List *outlist(new Particle_List());
  outlist->resize(inlist->size());
  for (size_t i(0);i<outlist->size();++i) {
    (*outlist)[i] = new Particle(*(*inlist)[i]);
    Vec4D p((*outlist)[i]->Momentum());
    zrot.Rotate(p);
    (*outlist)[i]->SetMomentum(p);
  } 
  p_ana->AddParticleList(m_outlist,outlist);
}

void Rotator::Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout)
{
  CreateParticleList();
}

Analysis_Object *Rotator::GetCopy() const
{
  return new Rotator(m_inlist,m_reflist,m_outlist,m_flavs,m_items);
}


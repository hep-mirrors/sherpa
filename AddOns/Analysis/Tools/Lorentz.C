#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
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
    void Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount);
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
    void Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount);
    Analysis_Object *GetCopy() const;    
  };// end of class Rotator

  class RBooster : public Analysis_Object {
  private:
    std::string  m_inlist, m_ref, m_outlist;
  public:
    RBooster(const std::string &inlist,
	    const std::string &ref,
	    const std::string &outlist);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount);
    Analysis_Object *GetCopy() const;    
  };// end of class RBooster

  class RRotator : public Analysis_Object {
  private:
    std::string  m_inlist, m_ref, m_outlist;
  public:
    RRotator(const std::string &inlist,
	    const std::string &ref,
	    const std::string &outlist);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount);
    Analysis_Object *GetCopy() const;    
  };// end of class RRotator

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(Booster,"CMSBoost",
 	       Analysis_Object,Argument_Matrix);

void ATOOLS::Getter<Analysis_Object,Argument_Matrix,Booster>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"RefList list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Flavs   flav1 .. flavN\n"
     <<std::setw(width+7)<<" "<<"Items   item1 .. itemN\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *ATOOLS::Getter<Analysis_Object,Argument_Matrix,Booster>::
operator()(const Argument_Matrix &parameters) const
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
	flavs.push_back(Flavour((kf_code)abs(kf)));
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

DECLARE_GETTER(Rotator,"ZRotate",
 	       Analysis_Object,Argument_Matrix);

void ATOOLS::Getter<Analysis_Object,Argument_Matrix,Rotator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"RefList list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Flavs   flav1 .. flavN\n"
     <<std::setw(width+7)<<" "<<"Items   item1 .. itemN\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *ATOOLS::Getter<Analysis_Object,Argument_Matrix,Rotator>::
operator()(const Argument_Matrix &parameters) const
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
	flavs.push_back(Flavour((kf_code)abs(kf)));
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

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

Booster::Booster(const std::string &inlist,
		 const std::string &reflist,
		 const std::string &outlist,
		 const std::vector<Flavour> &flavs,
		 const std::vector<int> &items):
  m_inlist(inlist), m_reflist(reflist), m_outlist(outlist), 
  m_flavs(flavs), m_items(items) 
{
  m_name="Boost";
}

void Booster::CreateParticleList()
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  Particle_List *reflist(p_ana->GetParticleList(m_reflist));
  if (inlist==NULL || reflist==NULL) {
    msg_Error()<<METHOD<<"(): Missing lists: '"<<m_inlist
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

void Booster::Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount)
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
  m_flavs(flavs), m_items(items) 
{
  m_name="Rot";
}

void Rotator::CreateParticleList()
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  Particle_List *reflist(p_ana->GetParticleList(m_reflist));
  if (inlist==NULL || reflist==NULL) {
    msg_Error()<<METHOD<<"(): Missing lists: '"<<m_inlist
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

void Rotator::Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount)
{
  CreateParticleList();
}

Analysis_Object *Rotator::GetCopy() const
{
  return new Rotator(m_inlist,m_reflist,m_outlist,m_flavs,m_items);
}

DECLARE_GETTER(RBooster,"CMSRBoost",
 	       Analysis_Object,Argument_Matrix);

void ATOOLS::Getter<Analysis_Object,Argument_Matrix,RBooster>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"RefMom  tag\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *ATOOLS::Getter<Analysis_Object,Argument_Matrix,RBooster>::
operator()(const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Selected"), reflist("FinalState");
  std::vector<Flavour> flavs;
  std::vector<int> items;
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="RefMom" && cur.size()>1) reflist=cur[1];
  }
  items.resize(flavs.size(),0);
  return new RBooster(inlist,reflist,outlist);
}

DECLARE_GETTER(RRotator,"ZRRotate",
 	       Analysis_Object,Argument_Matrix);

void ATOOLS::Getter<Analysis_Object,Argument_Matrix,RRotator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"RefMom  tag\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *ATOOLS::Getter<Analysis_Object,Argument_Matrix,RRotator>::
operator()(const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Selected"), reflist("FinalState");
  std::vector<Flavour> flavs;
  std::vector<int> items;
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="RefMom" && cur.size()>1) reflist=cur[1];
  }
  items.resize(flavs.size(),0);
  return new RRotator(inlist,reflist,outlist);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

RBooster::RBooster(const std::string &inlist,
		 const std::string &ref,
		 const std::string &outlist):
  m_inlist(inlist), m_ref(ref), m_outlist(outlist) 
{
  m_name="RBoost";
}

void RBooster::CreateParticleList()
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  Blob_Data_Base *data((*p_ana)[m_ref]);
  if (data==NULL) THROW(fatal_error,"Reference momentum not found");
  Vec4D cms(data->Get<Vec4D>());
  Poincare cmsrboost(cms);
  Particle_List *outlist(new Particle_List());
  outlist->resize(inlist->size());
  for (size_t i(0);i<outlist->size();++i) {
    (*outlist)[i] = new Particle(*(*inlist)[i]);
    Vec4D p((*outlist)[i]->Momentum());
    cmsrboost.Boost(p);
    (*outlist)[i]->SetMomentum(p);
  } 
  p_ana->AddParticleList(m_outlist,outlist);
}

void RBooster::Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount)
{
  CreateParticleList();
}

Analysis_Object *RBooster::GetCopy() const
{
  return new RBooster(m_inlist,m_ref,m_outlist);
}

RRotator::RRotator(const std::string &inlist,
		 const std::string &ref,
		 const std::string &outlist):
  m_inlist(inlist), m_ref(ref), m_outlist(outlist) 
{
  m_name="RRot";
}

void RRotator::CreateParticleList()
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  Blob_Data_Base *data((*p_ana)[m_ref]);
  if (data==NULL) THROW(fatal_error,"Reference momentum not found");
  Vec4D zaxis(data->Get<Vec4D>());
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

void RRotator::Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount)
{
  CreateParticleList();
}

Analysis_Object *RRotator::GetCopy() const
{
  return new RRotator(m_inlist,m_ref,m_outlist);
}


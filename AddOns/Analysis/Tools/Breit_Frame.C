#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Breit_Frame : public Analysis_Object {
  private:
    std::string  m_inlist, m_outlist;
  public:
    Breit_Frame(const std::string &inlist,
	    const std::string &outlist);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount);
    Analysis_Object *GetCopy() const;    
  };// end of class Breit_Frame

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(Breit_Frame_Getter,"BreitFrame",
 	       Analysis_Object,Argument_Matrix);

void Breit_Frame_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object *
Breit_Frame_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Selected");
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
  }
  return new Breit_Frame(inlist,outlist);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

Breit_Frame::Breit_Frame(const std::string &inlist,
		       const std::string &outlist):
  m_inlist(inlist), m_outlist(outlist)
{
  m_name="Breit_Frame";
}

void Breit_Frame::Evaluate(const ATOOLS::Blob_List &bl,double weight,double ncount)
{
  Vec4D l, lp, pp;
  for (size_t i(0);i<bl.size();++i)
    if (bl[i]->Type()==btp::Beam) {
      Particle *p(bl[i]->InParticle(0));
      if (p->Flav().IsLepton()) l=p->Momentum();
      else pp=p->Momentum();
    }
  Blob *me(bl.FindFirst(btp::Signal_Process));
  for (int i(0);i<me->NOutP();++i)
    if (me->OutParticle(i)->Flav().IsLepton()) {
      lp=me->OutParticle(i)->Momentum();
      break;
    }
  Vec4D qq(l-lp);
  Poincare cms(pp+qq);
  double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), E(sqrt(Q2)/(2.0*x));
  Vec4D mp=Vec4D(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
  Vec4D mq=Vec4D(0.0,0.0,0.0,2.0*x*E);
  cms.Boost(pp);
  cms.Boost(qq);
  Poincare zrot(pp,-Vec4D::ZVEC);
  zrot.Rotate(pp);
  zrot.Rotate(qq);
  Poincare breit(mp+mq);
  breit.BoostBack(pp);
  breit.BoostBack(qq);
  if (!IsEqual(pp,mp,1.0e-3) || !IsEqual(qq,mq,1.0e-3))
    msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg_Error()<<METHOD<<"(): Missing list: '"
	       <<m_inlist<<"'."<<std::endl;
  }
  Particle_List *outlist(new Particle_List());
  outlist->resize(inlist->size());
  for (size_t i(0);i<outlist->size();++i) {
    (*outlist)[i] = new Particle(*(*inlist)[i]);
    Vec4D p((*outlist)[i]->Momentum());
    cms.Boost(p);
    zrot.Rotate(p);
    breit.BoostBack(p);
    (*outlist)[i]->SetMomentum(p);
  } 
  p_ana->AddParticleList(m_outlist,outlist);
}

Analysis_Object *Breit_Frame::GetCopy() const
{
  return new Breit_Frame(m_inlist,m_outlist);
}

#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Phys/Ordering.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <algorithm>

namespace ANALYSIS {

  class Tau_Reconstruction: public Trigger_Base {  
  private:

    ATOOLS::Order_Base *p_order;

  public:

    Tau_Reconstruction(const std::string &inlist,const std::string &outlist,
		       Primitive_Analysis *const ana);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double weight, double ncount);
    Analysis_Object *GetCopy() const;    
    Analysis_Object &operator+=(const Analysis_Object & ob);

    void Output(const std::string & pname);

  };// end of class Tau_Reconstruction

}// end of namespace ANALYSIS

using namespace ANALYSIS;
using namespace ATOOLS;

template <class Class>
Analysis_Object *GetTrigger(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters[0].size()<2) return NULL;
  return new Class(parameters[0][0],parameters[0][1],parameters());
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Analysis_Object *							\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetTrigger<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"inlist outlist"; }

#define DEFINE_TRIGGER_GETTER(CLASS,NAME,TAG)				\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);		\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_TRIGGER_GETTER(Tau_Reconstruction,Tau_Reconstruction_Getter,"TauReco")

Tau_Reconstruction::Tau_Reconstruction(const std::string &inlist,
				       const std::string &outlist,
				       Primitive_Analysis *const ana):
  Trigger_Base(inlist,outlist)
{
  p_order = Order_Getter::GetObject("PT_UP","");
  if (p_order==NULL) THROW(fatal_error,"Invalid ordering mode.");
  ana->AddData(m_outlist+"_x1",new Blob_Data<double>(0.0));
  ana->AddData(m_outlist+"_x2",new Blob_Data<double>(0.0));
  ana->AddData(m_outlist+"_mtautau",new Blob_Data<double>(0.0));  
}

void Tau_Reconstruction::Evaluate(const ATOOLS::Particle_List &inlist,
				  ATOOLS::Particle_List &outlist,
				  double weight, double ncount)
{
  Particle *l1(NULL), *l2(NULL), *etmiss(NULL);
  Particle_List mylist(inlist);
  std::sort(mylist.begin(),mylist.end(),*p_order);
  for (Particle_List::const_iterator pit=mylist.begin();
       pit!=mylist.end();++pit) {
    if ((*pit)->Flav().IsLepton()) {
      if (l2!=NULL) {
        msg_Error()<<METHOD<<"(): More than two hard leptons in '"
                   <<m_inlist<<"'"<<std::endl;
      }
      if (l1!=NULL) l2=*pit;
      else l1=*pit;
    }
    else if ((*pit)->Flav().Kfcode()==kf_none) {
      if (etmiss!=NULL) {
        msg_Error()<<METHOD<<"(): More than one missing et particle in '"
                   <<m_inlist<<"'"<<std::endl;
      }
      etmiss=*pit;
    }
    if (l1 && l2 && etmiss) break;
  }
  if (l1==NULL || l2==NULL || etmiss==NULL) return;
  Vec3D ptl1(l1->Momentum().Perp());
  Vec3D ptl2(l2->Momentum().Perp());
  Vec3D ptet(etmiss->Momentum().Perp());
  double num(cross(ptl1,ptl2)[3]);
  double x1(num/(num+cross(ptet,ptl2)[3])); 
  double x2(num/(num+cross(ptl1,ptet)[3]));
  double mtautau((l1->Momentum()+l2->Momentum()).Mass()/sqrt(x1*x2));
  // put into ana
  p_ana->AddData(m_outlist+"_x1",new Blob_Data<double>(x1));
  p_ana->AddData(m_outlist+"_x2",new Blob_Data<double>(x2));
  p_ana->AddData(m_outlist+"_mtautau",new Blob_Data<double>(mtautau));  
}

Analysis_Object &Tau_Reconstruction::operator+=(const Analysis_Object & ob)
{
  return *this;
}

void Tau_Reconstruction::Output(const std::string & pname) 
{
}

Analysis_Object *Tau_Reconstruction::GetCopy() const 
{
  return new Tau_Reconstruction(m_inlist,m_outlist,p_ana);
}

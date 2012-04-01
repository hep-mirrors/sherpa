#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Four_Particle_Calculator_Base: public Analysis_Object {  
  protected:

    std::string m_inlist, m_outlist;

    double m_xmin, m_xmax;

    ATOOLS::Flavour m_flav[4];
    size_t          m_item[4];

  public:

    Four_Particle_Calculator_Base
    (const ATOOLS::Flavour flav[4],const size_t item[4],
     const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Blob_List &bl,
		  double weight,double ncount);
    
    virtual bool Calculate(const Particle *p1,const Particle *p2,
			   const Particle *p3,const Particle *p4) const = 0;

  };// end of class Four_Particle_Calculator_Base

  class Intermediate_Mass_Fitter: public Four_Particle_Calculator_Base {  
  public:

    Intermediate_Mass_Fitter
    (const ATOOLS::Flavour flav[4],const size_t item[4],
     const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    bool Calculate(const Particle *p1,const Particle *p2,
		   const Particle *p3,const Particle *p4) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Intermediate_Mass_Fitter

  class Delta_Phi_Fitter: public Four_Particle_Calculator_Base {  
  public:

    Delta_Phi_Fitter
    (const ATOOLS::Flavour flav[4],const size_t item[4],
     const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    bool Calculate(const Particle *p1,const Particle *p2,
		   const Particle *p3,const Particle *p4) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Delta_Phi_Fitter

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetFourParticleCalculator(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<12) return NULL;
    size_t item[4];
    ATOOLS::Flavour flav[4];
    for (size_t i(0);i<4;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][2*i]);
      flav[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) flav[i]=flav[i].Bar();
      item[i]=ATOOLS::ToType<size_t>(parameters[0][2*i+1]);
    }
    return new Class(flav,item,
		     ATOOLS::ToType<double>(parameters[0][8]),
		     ATOOLS::ToType<double>(parameters[0][9]),
		     parameters[0][10],parameters[0][11]);
  }
  return NULL;
}									

#define DEFINE_FOUR_CALCULATOR_GETTER_METHOD(CLASS,NAME)		\
  Analysis_Object *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetFourParticleCalculator<CLASS>(parameters); }

#define DEFINE_FOUR_CALCULATOR_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav1 item1 ... flav4 item4 min max inlist resulttag"; }

#define DEFINE_FOUR_CALCULATOR_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_FOUR_CALCULATOR_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_FOUR_CALCULATOR_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Four_Particle_Calculator_Base::Four_Particle_Calculator_Base    
(const ATOOLS::Flavour flav[4],const size_t item[4],
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  m_inlist(inlist), m_outlist(outlist)
{
  for (size_t i(0);i<4;++i) {
    m_flav[i]=flav[i];
    m_item[i]=item[i];
  }
  m_xmin=min;
  m_xmax=max;
}

void Four_Particle_Calculator_Base::Evaluate
(const ATOOLS::Blob_List &bl,double weight,double ncount)
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg_Error()<<METHOD<<"(): List '"<<m_inlist
		       <<"' not found."<<std::endl;
    return;
  }
  int no(-1);
  size_t pos[4]={std::string::npos,std::string::npos,
		 std::string::npos,std::string::npos};
  for (size_t k(0);k<4;++k) {
    no=-1;
    for (size_t i(0);i<inlist->size();++i) {
      if ((*inlist)[i]->Flav()==m_flav[k] || 
	  m_flav[k].Kfcode()==kf_none) {
	++no;
	if (no==(int)m_item[k]) {
	  pos[k]=i;
	  break;
	}
      }
    }
  }
  for (size_t k(0);k<4;++k) if (pos[k]==std::string::npos) return;
  Calculate((*inlist)[pos[0]],(*inlist)[pos[1]],
	    (*inlist)[pos[2]],(*inlist)[pos[3]]);
}

DECLARE_GETTER(IMF_Getter,"BWMassFit",
	       Analysis_Object,Argument_Matrix);
Analysis_Object *					
IMF_Getter::operator()(const Argument_Matrix &parameters) const	       
{ return GetFourParticleCalculator<Intermediate_Mass_Fitter>(parameters); }
void IMF_Getter::PrintInfo(std::ostream &str,const size_t width) const       
{ str<<"flav1 item1 ... flav4 item4 mass width inlist resulttag"; }
  
Intermediate_Mass_Fitter::Intermediate_Mass_Fitter
(const ATOOLS::Flavour flav[4],const size_t item[4],
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Four_Particle_Calculator_Base(flav,item,min,max,inlist,outlist) {}

bool Intermediate_Mass_Fitter::
Calculate(const Particle *p1,const Particle *p2,
	  const Particle *p3,const Particle *p4) const
{
  // 1st combination
  double mass12((p1->Momentum()+p2->Momentum()).Abs2());
  double mass34((p3->Momentum()+p4->Momentum()).Abs2());
  double pr12(1.0/(sqr(mass12-m_xmin*m_xmin)+sqr(m_xmin*m_xmax)));
  double pr34(1.0/(sqr(mass34-m_xmin*m_xmin)+sqr(m_xmin*m_xmax)));
  double pr1(pr12*pr34);
  // 2nd combination
  double mass14((p1->Momentum()+p4->Momentum()).Abs2());
  double mass32((p3->Momentum()+p2->Momentum()).Abs2());
  double pr14(1.0/(sqr(mass14-m_xmin*m_xmin)+sqr(m_xmin*m_xmax)));
  double pr32(1.0/(sqr(mass32-m_xmin*m_xmin)+sqr(m_xmin*m_xmax)));
  double pr2(pr14*pr32);
  msg_Debugging()<<"m12 = "<<sqrt(mass12)<<",m14 = "<<sqrt(mass14)
		 <<",m34 = "<<sqrt(mass34)<<",m32 = "<<sqrt(mass32)
		 <<" -> add '"<<m_outlist<<"' = "<<log10(pr1/pr2)<<"\n";
  p_ana->AddData(m_outlist,new Blob_Data<double>(log10(pr1/pr2)));
  return true;
}

Analysis_Object *Intermediate_Mass_Fitter::GetCopy() const
{
  return new Intermediate_Mass_Fitter(m_flav,m_item,m_xmin,m_xmax,
				      m_inlist,m_outlist);
}

DECLARE_GETTER(DPF_Getter,"CosDPhiFit",
	       Analysis_Object,Argument_Matrix);
Analysis_Object *					
DPF_Getter::operator()(const Argument_Matrix &parameters) const	       
{ return GetFourParticleCalculator<Delta_Phi_Fitter>(parameters); }
void DPF_Getter::PrintInfo(std::ostream &str,const size_t width) const       
{ str<<"flav1 item1 ... flav4 item4 dphi1 dphi2 inlist resulttag"; }
  
Delta_Phi_Fitter::Delta_Phi_Fitter
(const ATOOLS::Flavour flav[4],const size_t item[4],
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Four_Particle_Calculator_Base(flav,item,min,max,inlist,outlist) {}

bool Delta_Phi_Fitter::
Calculate(const Particle *p1,const Particle *p2,
	  const Particle *p3,const Particle *p4) const
{
  // 1st combination
  double dphi12(p1->Momentum().DPhi(p2->Momentum()));
  double dphi34(p3->Momentum().DPhi(p4->Momentum()));
  double pr1(cos(dphi12-m_xmin)*cos(dphi34-m_xmax));
  // 2nd combination
  double dphi14(p1->Momentum().DPhi(p4->Momentum()));
  double dphi32(p3->Momentum().DPhi(p2->Momentum()));
  double pr2(cos(dphi14-m_xmin)*cos(dphi32-m_xmax));
  msg_Debugging()<<"dphi12 = "<<dphi12<<",dphi14 = "<<dphi14
		 <<",dphi34 = "<<dphi34<<",dphi32 = "<<dphi32
		 <<" -> add '"<<m_outlist<<"' = "<<log10(pr1/pr2)<<"\n";
  p_ana->AddData(m_outlist,new Blob_Data<double>(log10(pr1/pr2)));
  return true;
}

Analysis_Object *Delta_Phi_Fitter::GetCopy() const
{
  return new Delta_Phi_Fitter(m_flav,m_item,m_xmin,m_xmax,
			      m_inlist,m_outlist);
}


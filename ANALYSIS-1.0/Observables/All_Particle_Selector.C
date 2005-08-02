#include "All_Particle_Selector.H"

#include "MyStrStream.H"
#include <iomanip>

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const 
GetParticleSelector(const String_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    return new Class(ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     parameters[0][2],parameters[0][3]);
  }
  if (parameters.size()<4) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Analysed", outlist="Selected";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
  }
  return new Class(min,max,inlist,outlist);
}									

#define DEFINE_SELECTOR_GETTER_METHOD(CLASS,NAME)		\
  Primitive_Observable_Base *				\
  NAME::operator()(const String_Matrix &parameters) const	\
  { return GetParticleSelector<CLASS>(parameters); }

#define DEFINE_SELECTOR_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max inlist outlist"; }

#define DEFINE_SELECTOR_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_SELECTOR_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_SELECTOR_PRINT_METHOD(NAME)

#include "Primitive_Analysis.H"

template <class Class>
Primitive_Observable_Base *const 
GetParticleDSelector(const String_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][3]);
    ATOOLS::Flavour flav((ATOOLS::kf::code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     flav,ATOOLS::ToType<int>(parameters[0][3]),
		     parameters[0][5],parameters[0][4],parameters[0][6]);
  }
  if (parameters.size()<7) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Analysed", reflist="Jets", outlist="Selected";
  size_t item=0;
  ATOOLS::Flavour flav(ATOOLS::kf::jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour(ATOOLS::kf::code(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
  }
  return new Class(min,max,flav,item,reflist,inlist,outlist);
}									

#define DEFINE_SELECTOR_D_GETTER_METHOD(CLASS,NAME)		\
  Primitive_Observable_Base *				\
  NAME::operator()(const String_Matrix &parameters) const	\
  { return GetParticleDSelector<CLASS>(parameters); }

#define DEFINE_SELECTOR_D_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max flav item inlist reflist outlist"; }

#define DEFINE_SELECTOR_D_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_SELECTOR_D_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_SELECTOR_D_PRINT_METHOD(NAME)

#include "Primitive_Analysis.H"

DEFINE_SELECTOR_GETTER(PT_Selector,PT_Selector_Getter,"PTSel")

PT_Selector::
PT_Selector(const double min,const double max,
	    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<PT<"+
	    ATOOLS::ToString(max)+inlist)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void PT_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  for (size_t i=0;i<particlelist.size();++i) {
    double pt=particlelist[i]->Momentum().PPerp();
    if (pt>=m_xmin && pt<=m_xmax) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
  }
}

Primitive_Observable_Base *PT_Selector::Copy() const
{
  return new PT_Selector(m_xmin,m_xmax,m_listname,m_outlist);
}

void PT_Selector::EndEvaluation(double scale)
{
}

DEFINE_SELECTOR_GETTER(ET_Selector,ET_Selector_Getter,"ETSel")

ET_Selector::
ET_Selector(const double min,const double max,
	    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<ET<"+
	    ATOOLS::ToString(max)+inlist)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void ET_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  for (size_t i=0;i<particlelist.size();++i) {
    double et=particlelist[i]->Momentum().EPerp();
    if (et>=m_xmin && et<=m_xmax) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
  }
}

Primitive_Observable_Base *ET_Selector::Copy() const
{
  return new ET_Selector(m_xmin,m_xmax,m_listname,m_outlist);
}

void ET_Selector::EndEvaluation(double scale)
{
}

DEFINE_SELECTOR_GETTER(Eta_Selector,Eta_Selector_Getter,"EtaSel")

Eta_Selector::
Eta_Selector(const double min,const double max,
	    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<Eta<"+
	    ATOOLS::ToString(max)+inlist)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void Eta_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  for (size_t i=0;i<particlelist.size();++i) {
    double eta=particlelist[i]->Momentum().Eta();
    if (eta>=m_xmin && eta<=m_xmax) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
  }
}

Primitive_Observable_Base *Eta_Selector::Copy() const
{
  return new Eta_Selector(m_xmin,m_xmax,m_listname,m_outlist);
}

void Eta_Selector::EndEvaluation(double scale)
{
}

DEFINE_SELECTOR_GETTER(Y_Selector,Y_Selector_Getter,"YSel")

Y_Selector::
Y_Selector(const double min,const double max,
	    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<Y<"+
	    ATOOLS::ToString(max)+inlist)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void Y_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  for (size_t i=0;i<particlelist.size();++i) {
    double y=particlelist[i]->Momentum().Y();
    if (y>=m_xmin && y<=m_xmax) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
  }
}

Primitive_Observable_Base *Y_Selector::Copy() const
{
  return new Eta_Selector(m_xmin,m_xmax,m_listname,m_outlist);
}

void Y_Selector::EndEvaluation(double scale)
{
}

DEFINE_SELECTOR_GETTER(Phi_Selector,Phi_Selector_Getter,"PhiSel")

Phi_Selector::
Phi_Selector(const double min,const double max,
	    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<Phi<"+
	    ATOOLS::ToString(max)+inlist)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void Phi_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  for (size_t i=0;i<particlelist.size();++i) {
    double phi=particlelist[i]->Momentum().Phi()/M_PI*180.0;
    if (phi>=m_xmin && phi<=m_xmax) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
  }
}

Primitive_Observable_Base *Phi_Selector::Copy() const
{
  return new Phi_Selector(m_xmin,m_xmax,m_listname,m_outlist);
}

void Phi_Selector::EndEvaluation(double scale)
{
}

DEFINE_SELECTOR_D_GETTER(DPhi_Selector,DPhi_Selector_Getter,"DPhiSel")

DPhi_Selector::
DPhi_Selector(const double min,const double max,
	      const ATOOLS::Flavour flav,
	      const size_t item,const std::string &reflist,
	      const std::string &inlist,const std::string &outlist):
  m_reflist(reflist),
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<DPhi<"+
	    ATOOLS::ToString(max)+inlist),
  m_item(item),
  m_flavour(flav)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void DPhi_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<reflist->size();++i) 
    if ((*reflist)[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  for (size_t i=0;i<particlelist.size();++i) {
    double phi=
      ATOOLS::dabs(particlelist[i]->Momentum().
		   DPhi((*reflist)[pos]->Momentum())/M_PI*180.0);
    if (phi>=m_xmin && phi<=m_xmax) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
  }
}

Primitive_Observable_Base *DPhi_Selector::Copy() const
{
  return new DPhi_Selector(m_xmin,m_xmax,m_flavour,m_item,
			   m_reflist,m_listname,m_outlist);
}

void DPhi_Selector::EndEvaluation(double scale)
{
}

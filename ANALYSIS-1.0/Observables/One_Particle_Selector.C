#include "One_Particle_Selector.H"

#include "MyStrStream.H"
#include <iomanip>

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const 
GetOneParticleSelector(const String_Matrix &parameters) 
{				
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((ATOOLS::kf::code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     parameters[0][5],parameters[0][6]);
  }
  if (parameters.size()<7) return NULL;
  double min=30.0, max=70.0; 
  std::string inlist="Jets", outlist="LeadJets";
  size_t item=0;
  int mode=0;
  ATOOLS::Flavour flav(ATOOLS::kf::jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Mode") mode=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour(ATOOLS::kf::code(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
  }
  return new Class(flav,item,mode,min,max,inlist,outlist);
}									

#define DEFINE_ONE_SELECTOR_GETTER_METHOD(CLASS,NAME)	\
  Primitive_Observable_Base *				\
  NAME::operator()(const String_Matrix &parameters) const	\
  { return GetOneParticleSelector<CLASS>(parameters); }

#define DEFINE_ONE_SELECTOR_PRINT_METHOD(NAME)			\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav item mode min max inlist outlist"; }

#define DEFINE_ONE_SELECTOR_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_ONE_SELECTOR_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_ONE_SELECTOR_PRINT_METHOD(NAME)

template <class Class>
Primitive_Observable_Base *const 
GetOneParticleDeltaSelector(const String_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<9) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((ATOOLS::kf::code)abs(kf));
    if (kf<0) flav=flav.Bar();
    kf=ATOOLS::ToType<int>(parameters[0][2]);
    ATOOLS::Flavour refflav((ATOOLS::kf::code)abs(kf));
    if (kf<0) refflav=refflav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     refflav,ATOOLS::ToType<size_t>(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][5]),
		     parameters[0][6],parameters[0][7],parameters[0][8]);
  }
  if (parameters.size()<9) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", reflist="Jets", outlist="LeadJets";
  size_t item=0, refitem=1;
  ATOOLS::Flavour flav(ATOOLS::kf::jet), refflav(ATOOLS::kf::jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item1") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Item2") refitem=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav1") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour(ATOOLS::kf::code(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
    else if (parameters[i][0]=="Flav2") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      refflav=ATOOLS::Flavour(ATOOLS::kf::code(abs(kf)));
      if (kf<0) refflav=refflav.Bar();
    }
  }
  return new Class(flav,item,refflav,refitem,min,max,inlist,reflist,outlist);
}									

#define DEFINE_ONE_SELECTOR_DELTA_GETTER_METHOD(CLASS,NAME)		\
  Primitive_Observable_Base *					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetOneParticleDeltaSelector<CLASS>(parameters); }

#define DEFINE_ONE_SELECTOR_DELTA_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav1 item1 flav2 item2 min max inlist reflist outlist"; }

#define DEFINE_ONE_SELECTOR_DELTA_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_ONE_SELECTOR_DELTA_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_ONE_SELECTOR_DELTA_PRINT_METHOD(NAME)

#include "Primitive_Analysis.H"

DEFINE_ONE_SELECTOR_GETTER(One_PT_Selector,
			   One_PT_Selector_Getter,"OnePTSel")

One_PT_Selector::
One_PT_Selector(const ATOOLS::Flavour flav,
		const size_t item,const int mode,const double min,const double max,
		const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_PT<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav), m_item(item), m_mode(mode)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_PT_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<particlelist.size();++i) 
    if (particlelist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  double pt=particlelist[pos]->Momentum().PPerp();
  if (pt<m_xmin || pt>m_xmax) return;
  
  if (m_mode==0) { 
    outlist->resize(particlelist.size());
  for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
  }
  else {
    int diff_flavour=0;
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour) ++diff_flavour;
    }

    int size = diff_flavour+1; 
    outlist->resize(size);
    
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour || i==pos) {
	(*outlist)[size-1] = new ATOOLS::Particle(*particlelist[i]);
	--size;
      }
    }
  }
}

Primitive_Observable_Base *One_PT_Selector::Copy() const
{
  return new One_PT_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,m_listname,m_outlist);
}

void One_PT_Selector::EndEvaluation(double scale)
{
}

DEFINE_ONE_SELECTOR_GETTER(One_ET_Selector,
			   One_ET_Selector_Getter,"OneETSel")

One_ET_Selector::
One_ET_Selector(const ATOOLS::Flavour flav,
		const size_t item,const int mode,const double min,const double max,
		const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_ET<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav), m_item(item), m_mode(mode)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_ET_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<particlelist.size();++i) 
    if (particlelist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  double et=particlelist[pos]->Momentum().EPerp();
  if (et<m_xmin || et>m_xmax) return;
  if (m_mode==0) { 
    outlist->resize(particlelist.size());
    for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
  }
  else {
    int diff_flavour=0;
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour) ++diff_flavour;
    }
    
    int size = diff_flavour+1; 
    outlist->resize(size);
    
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour || i==pos) {
	(*outlist)[size-1] = new ATOOLS::Particle(*particlelist[i]);
	--size;
      }
    }
  }
}

Primitive_Observable_Base *One_ET_Selector::Copy() const
{
  return new One_ET_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,m_listname,m_outlist);
}

void One_ET_Selector::EndEvaluation(double scale)
{
}

DEFINE_ONE_SELECTOR_GETTER(One_Eta_Selector,
			   One_Eta_Selector_Getter,"OneEtaSel")

One_Eta_Selector::
One_Eta_Selector(const ATOOLS::Flavour flav,
		const size_t item,const int mode,const double min,const double max,
		const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_Eta<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav), m_item(item), m_mode(mode)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_Eta_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<particlelist.size();++i) 
    if (particlelist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  double eta=particlelist[pos]->Momentum().Eta();
  if (eta<m_xmin || eta>m_xmax) return;
  
  if (m_mode==0) { 
    outlist->resize(particlelist.size());
    for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
  }
  else {
    int diff_flavour=0;
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour) ++diff_flavour;
    }
    
    int size = diff_flavour+1; 
    outlist->resize(size);
    
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour || i==pos) {
	(*outlist)[size-1] = new ATOOLS::Particle(*particlelist[i]);
	--size;
      }
    }
  }
}

Primitive_Observable_Base *One_Eta_Selector::Copy() const
{
  return new One_Eta_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,m_listname,m_outlist);
}

void One_Eta_Selector::EndEvaluation(double scale)
{
}

DEFINE_ONE_SELECTOR_GETTER(One_Y_Selector,
			   One_Y_Selector_Getter,"OneYSel")

One_Y_Selector::
One_Y_Selector(const ATOOLS::Flavour flav,
		const size_t item,const int mode,const double min,const double max,
		const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_Y<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav), m_item(item), m_mode(mode)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_Y_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<particlelist.size();++i) 
    if (particlelist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  double et=particlelist[pos]->Momentum().Y();
  if (et<m_xmin || et>m_xmax) return;
  if (m_mode==0) { 
    outlist->resize(particlelist.size());
    for (size_t i=0;i<particlelist.size();++i) 
      (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
  }
  else {
    int diff_flavour=0;
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour) ++diff_flavour;
    }
    
    int size = diff_flavour+1; 
    outlist->resize(size);
    
    for (size_t i=0;i<particlelist.size();++i) {
      if (particlelist[i]->Flav()!=m_flavour || i==pos) {
	(*outlist)[size-1] = new ATOOLS::Particle(*particlelist[i]);
	--size;
      }
    }
  }  
}

Primitive_Observable_Base *One_Y_Selector::Copy() const
{
  return new One_Y_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,m_listname,m_outlist);
}

void One_Y_Selector::EndEvaluation(double scale)
{
}

DEFINE_ONE_SELECTOR_DELTA_GETTER(One_DPhi_Selector,
				 One_DPhi_Selector_Getter,"OneDPhiSel")

One_DPhi_Selector::
One_DPhi_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  m_reflist(reflist),
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_DPhi<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav),
  m_refflavour(refflav),
  m_item(item),
  m_refitem(refitem)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_DPhi_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1, refno=-1; 
  size_t pos=std::string::npos, refpos=std::string::npos;
  for (size_t i=0;i<reflist->size();++i) {
    if ((*reflist)[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	if (refpos!=std::string::npos) break;
      }
    }
    if ((*reflist)[i]->Flav()==m_refflavour || 
	m_refflavour.Kfcode()==ATOOLS::kf::none) {
      ++refno;
      if (refno==(int)m_refitem) {
	refpos=i;
	if (pos!=std::string::npos) break;
      }
    }
  }
  if (pos==std::string::npos || refpos==std::string::npos) return;
  double dphi=
    ATOOLS::dabs((*reflist)[pos]->Momentum().
		 DPhi((*reflist)[refpos]->Momentum())/M_PI*180.0);
  if (dphi<m_xmin || dphi>m_xmax) return;
  outlist->resize(particlelist.size());
  for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
}

Primitive_Observable_Base *One_DPhi_Selector::Copy() const
{
  return new One_DPhi_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_listname,m_reflist,m_outlist);
}

void One_DPhi_Selector::EndEvaluation(double scale)
{
}

DEFINE_ONE_SELECTOR_DELTA_GETTER(One_ETFrac_Selector,
				 One_ETFrac_Selector_Getter,"OneETFracSel")

One_ETFrac_Selector::
One_ETFrac_Selector(const ATOOLS::Flavour flav,const size_t item,
		   const ATOOLS::Flavour refflav,const size_t refitem,
		   const double min,const double max,
		   const std::string &inlist,const std::string &reflist,
		   const std::string &outlist):
  m_reflist(reflist),
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_ETFrac<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav),
  m_refflavour(refflav),
  m_item(item),
  m_refitem(refitem)
{
  m_splitt_flag = false;
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_ETFrac_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1, refno=-1; 
  size_t pos=std::string::npos, refpos=std::string::npos;
  for (size_t i=0;i<reflist->size();++i) {
    if ((*reflist)[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	if (refpos!=std::string::npos) break;
      }
    }
    if ((*reflist)[i]->Flav()==m_refflavour || 
	m_refflavour.Kfcode()==ATOOLS::kf::none) {
      ++refno;
      if (refno==(int)m_refitem) {
	refpos=i;
	if (pos!=std::string::npos) break;
      }
    }
  }
  if (pos==std::string::npos || refpos==std::string::npos) return;
  double efrac=(*reflist)[pos]->Momentum().EPerp()/
    (*reflist)[refpos]->Momentum().EPerp();
  if (efrac<m_xmin || efrac>m_xmax) return;
  outlist->resize(particlelist.size());
  for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
}

Primitive_Observable_Base *One_ETFrac_Selector::Copy() const
{
  return new One_ETFrac_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				 m_xmin,m_xmax,m_listname,m_reflist,m_outlist);
}

void One_ETFrac_Selector::EndEvaluation(double scale)
{
}


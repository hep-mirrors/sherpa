#include "One_Particle_Extractor.H"

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
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     parameters[0][4],parameters[0][5],parameters[0][6]);
  }
  if (parameters.size()<6) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", reflist="Jets", outlist="LeadJets";
  size_t item=0;
  ATOOLS::Flavour flav(ATOOLS::kf::jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour(ATOOLS::kf::code(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
  }
  return new Class(flav,item,min,max,inlist,reflist,outlist);
}									

#define DEFINE_ONE_EXTRACTOR_GETTER_METHOD(CLASS,NAME)	\
  Primitive_Observable_Base *const				\
  NAME::operator()(const String_Matrix &parameters) const	\
  { return GetOneParticleSelector<CLASS>(parameters); }

#define DEFINE_ONE_EXTRACTOR_PRINT_METHOD(NAME)			\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav item min max inlist outlist"; }

#define DEFINE_ONE_EXTRACTOR_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_ONE_EXTRACTOR_GETTER_METHOD(CLASS,NAME);			\
  DEFINE_ONE_EXTRACTOR_PRINT_METHOD(NAME)

#include "Primitive_Analysis.H"

DEFINE_ONE_EXTRACTOR_GETTER(One_PT_Extractor,
			    One_PT_Extractor_Getter,"OnePTExt");

One_PT_Extractor::
One_PT_Extractor(const ATOOLS::Flavour flav,
		     const size_t item,const double min,const double max,
		     const std::string &inlist,const std::string &reflist,
		     const std::string &outlist):
  m_reflist(reflist),
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(min)+"<One_PT<"+
	    ATOOLS::ToString(max)+inlist),
  m_flavour(flav),
  m_item(item)
{
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void One_PT_Extractor::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List(particlelist.size());
  for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
  p_ana->AddParticleList(m_outlist,outlist);
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1; 
  for (size_t i=0;i<reflist->size();++i) 
    if ((*reflist)[i]->Flav()==m_flavour || m_flavour.Kfcode()==ATOOLS::kf::none) {
      ++no;
      if (no==(int)m_item) {
	double pt=(*reflist)[i]->Momentum().PPerp();
	if (pt>m_xmin && pt<m_xmax) {
	  for (ATOOLS::Particle_List::iterator pit=outlist->begin();
	       pit!=outlist->end();++pit) {
	    if (*pit==(*reflist)[i]) outlist->erase(pit);
	    return;
	  }
	  msg_Tracking()<<"One_PT_Extractor::Evaluate(..): "
			<<"Cannot extract particle ("<<(*reflist)[i]<<") from '"
			<<m_listname<<"'"<<std::endl;
	}
	break;
      }
    }
}

Primitive_Observable_Base *One_PT_Extractor::Copy() const
{
  return new One_PT_Extractor(m_flavour,m_item,m_xmin,m_xmax,
				  m_listname,m_reflist,m_outlist);
}

void One_PT_Extractor::EndEvaluation(double scale)
{
}


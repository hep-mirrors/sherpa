#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

class Field_Stuart_Cone_Finder: public Trigger_Base {
private:
  
  double m_pt2min, m_etamin, m_etamax, m_deltar, m_deltar2;

  int    m_rscheme;

public:

  // constructor
  Field_Stuart_Cone_Finder(const double &ptmin,const double &etamin,
			   const double &etamax,const double &deltar,
			   const std::string &inlist,
			   const std::string &outlist,const int rscheme):
    Trigger_Base(inlist,outlist),
    m_pt2min(ptmin*ptmin), m_etamin(etamin), m_etamax(etamax), 
    m_deltar(deltar), m_deltar2(deltar*deltar), m_rscheme(rscheme) {}

  // member functions
  Analysis_Object *GetCopy() const 
  {
    return new Field_Stuart_Cone_Finder
      (sqrt(m_pt2min),m_etamin,m_etamax,sqrt(m_deltar2),
       m_inlist,m_outlist,m_rscheme);
  }

  double DeltaR2(const ATOOLS::Particle *p1,const ATOOLS::Particle *p2)
  {
    const ATOOLS::Vec4D &m1=p1->Momentum(), &m2=p2->Momentum();
    return ATOOLS::sqr(m1.DEta(m2))+ATOOLS::sqr(m1.DPhi(m2));
  }

  void Combine(ATOOLS::Particle *const p1,ATOOLS::Particle *const p2)
  {
    switch (m_rscheme) {
    case 0:
      // E scheme
      p1->SetMomentum(p1->Momentum()+p2->Momentum());
      break;
    case 1: {
      // scalar P_T scheme
      ATOOLS::Vec4D jet(p1->Momentum()+p2->Momentum());
      jet*=1.0/jet.PSpat();
      jet[0]=1.0;
      p1->SetMomentum((p1->Momentum().PPerp()+p2->Momentum().PPerp())/
		      ATOOLS::dabs(jet.SinTheta())*jet);
      break;
    }
    default:
      THROW(critical_error,"Recombination scheme invalid");
    }
  }

  void Evaluate(const ATOOLS::Particle_List &particlelist,
		ATOOLS::Particle_List &outlist,
		double value,double ncount)
  {
    ATOOLS::Particle_List *jetlist(&outlist), addlist;
    // select particles
    double etaplus=m_etamin-m_deltar, etaminus=m_etamax+m_deltar;
    for (size_t i=0;i<particlelist.size();++i) {
      double eta=particlelist[i]->Momentum().Eta();
      if (particlelist[i]->Momentum().PPerp2()>=m_pt2min) {
	if (eta>=m_etamin && eta<=m_etamax) {
	  jetlist->push_back(new ATOOLS::Particle(*particlelist[i]));
	  jetlist->back()->SetFlav(kf_jet);
	}
	else if (eta>=etaplus && eta<=etaminus) {
	  addlist.push_back(new ATOOLS::Particle(*particlelist[i]));
	}
      }
    }
    // sort according to pt
    std::stable_sort(jetlist->begin(),jetlist->end(),Order_PT());
    std::stable_sort(addlist.begin(),addlist.end(),Order_PT());
    // build jets
    for (size_t i=0;i<jetlist->size();++i) {
      if ((*jetlist)[i]==NULL) continue;
      for (size_t j=i+1;j<jetlist->size();++j) {
	if ((*jetlist)[j]!=NULL && 
	    DeltaR2((*jetlist)[i],(*jetlist)[j])<m_deltar2) {
	  Combine((*jetlist)[i],(*jetlist)[j]);
	  delete (*jetlist)[j];
	  (*jetlist)[j]=NULL;
	}
      }
      for (size_t j=0;j<addlist.size();++j) {
	if (addlist[j]!=NULL && 
	    DeltaR2((*jetlist)[i],addlist[j])<m_deltar2) {
	  Combine((*jetlist)[i],addlist[j]);
	  delete addlist[j];
	  addlist[j]=NULL;
	}
      }
    }
    // delete obsolete
    for (ATOOLS::Particle_List::iterator pit=jetlist->begin();
	 pit!=jetlist->end();) {
      if (*pit==NULL) {
	jetlist->erase(pit);
	pit=jetlist->begin();
      }
      else {
	++pit;
      }
    }
    for (ATOOLS::Particle_List::iterator pit=addlist.begin();
	 pit!=addlist.end();++pit) if (*pit!=NULL) delete *pit;
  }

};

DECLARE_GETTER(FS_Cone_Getter,"FSCone",
	       Analysis_Object,Argument_Matrix);	

Analysis_Object *FS_Cone_Getter::
operator()(const Argument_Matrix &parameters) const	
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<6) return NULL;
    return new 
      Field_Stuart_Cone_Finder(ATOOLS::ToType<double>(parameters[0][0]),
			       ATOOLS::ToType<double>(parameters[0][1]),
			       ATOOLS::ToType<double>(parameters[0][2]),
			       ATOOLS::ToType<double>(parameters[0][3]),
			       parameters[0][4],parameters[0][5],
			       parameters[0].size()>6?
			       ATOOLS::ToType<int>(parameters[0][6]):0);
  }
  else if (parameters.size()<6) return NULL;
  int recom=0;
  double ptmin=0.5, etamin=-1.0, etamax=1.0, deltar=0.7;
  std::string inlist="FinalState", outlist="LeadPartJets";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="PTMin") ptmin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="EtaMin") etamin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="EtaMax") etamax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="DeltaR") deltar=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Recom") recom=ATOOLS::ToType<int>(parameters[i][1]);
  }
  return new 
    Field_Stuart_Cone_Finder(ptmin,etamin,etamax,deltar,inlist,outlist,recom);
}									

void FS_Cone_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"ptmin etamin etamax deltar inlist outlist [recom]"; 
}


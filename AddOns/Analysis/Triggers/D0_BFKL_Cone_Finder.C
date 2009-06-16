#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

class Order_SET {
public:
  int operator()(ATOOLS::Particle * a, ATOOLS::Particle * b) {
    if (a==NULL || b==NULL) return 0;
    if (a->Momentum().EPerp() > b->Momentum().EPerp()) return 1;
    return 0;
  }
  int operator()(ATOOLS::Vec4D & a, ATOOLS::Vec4D & b) {
    if (a.EPerp() > b.EPerp()) return 1;
    return 0;
  }
};
  
class D0_BFKL_Cone_Finder: public Trigger_Base {
private:
  
  double m_etmin, m_deltarpre2, m_deltar2, m_deltarmin2;

public:

  // constructor
  D0_BFKL_Cone_Finder(const double &etmin,const double &deltarpre,
		      const double &deltar,const double &deltarmin,
		      const std::string &inlist,
		      const std::string &outlist):
    Trigger_Base(inlist,outlist),
    m_etmin(etmin), m_deltarpre2(deltarpre*deltarpre), 
    m_deltar2(deltar*deltar), m_deltarmin2(deltarmin*deltarmin) {}

  // member functions
  Analysis_Object *GetCopy() const 
  {
    return new D0_BFKL_Cone_Finder(m_etmin,sqrt(m_deltarpre2),sqrt(m_deltar2),
				   sqrt(m_deltarmin2),m_inlist,m_outlist);
  }

  double DeltaR2(const ATOOLS::Vec4D &m1,const ATOOLS::Vec4D &m2)
  {
    return sqr(m1.DEta(m2))+sqr(m1.DPhi(m2));
  }

  double DeltaR2(const ATOOLS::Particle *p1,const ATOOLS::Particle *p2)
  {
    return DeltaR2(p1->Momentum(),p2->Momentum());
  }

  ATOOLS::Vec4D ECenter(ATOOLS::Particle_List const &pl)
  {
    Vec4D sum;
    for (Particle_List::const_iterator pit(pl.begin());
	 pit!=pl.end();++pit) {
      Vec4D cp((*pit)->Momentum());
      sum+=cp[0]/cp.PSpat()*cp;
    }
    sum[0]=sum.PSpat();
    return sum;
  }

  ATOOLS::Vec4D ETCenter(ATOOLS::Particle_List const &pl)
  {
    Vec4D sum;
    for (Particle_List::const_iterator pit(pl.begin());
	 pit!=pl.end();++pit) {
      Vec4D cp((*pit)->Momentum());
      sum+=cp.EPerp()/cp.PSpat()*cp;
    }
    sum[0]=sum.PSpat();
    return sum;
  }

  bool Add(ATOOLS::Particle_List &pl,ATOOLS::Particle *const p)
  {
    for (Particle_List::iterator pit(pl.begin());pit!=pl.end();++pit)
      if (*pit==p) return false;
    pl.push_back(p);
    return true;
  }

  bool Remove(ATOOLS::Particle_List &pl,ATOOLS::Particle *const p)
  {
    for (Particle_List::iterator pit(pl.begin());pit!=pl.end();++pit)
      if (*pit==p) {
	pl.erase(pit);
	return true;
      }
    return false;
  }

  void Evaluate(const ATOOLS::Particle_List &inlist,
		ATOOLS::Particle_List &outlist,
		double weight, double ncount)
  {
    msg_Debugging()<<METHOD<<"() {\n";
    Particle_List towlist(inlist), stowlist(towlist);
    // sort in et
    std::stable_sort(towlist.begin(),towlist.end(),Order_ET());
    // construct preclusters
    for (Particle_List::iterator 
	   pit(towlist.begin());pit!=towlist.end();++pit) 
      msg_Debugging()<<"  tower "<<**pit<<"\n";
    msg_Debugging()<<"  preclustering {\n";
    while (towlist.size() && 
	   towlist.front()->Momentum().EPerp()>m_etmin) {
      Particle *seed(towlist.front());
      msg_Debugging()<<"    seed tower "<<*seed<<"\n";
      Particle_List tlist;
      for (Particle_List::iterator 
	     pit(towlist.begin());pit!=towlist.end();) {
	if ((*pit)->Momentum().EPerp()<m_etmin) break;
	if (DeltaR2(seed,*pit)>m_deltarpre2) {
	  ++pit;
	}
	else {
	  msg_Debugging()<<"    associated "<<**pit<<", dr="
			 <<sqrt(DeltaR2(seed,*pit))<<"\n";
	  tlist.push_back(*pit);
	  pit=towlist.erase(pit);
	}
      }
      outlist.push_back(new Particle(-1,kf_jet,ETCenter(tlist)));
      outlist.back()->SetInfo('P');
      msg_Debugging()<<"    precluster "<<*outlist.back()<<"\n";
    }
    msg_Debugging()<<"  }\n";
    if (outlist.empty()) return;
    // construct clusters
    std::map<Particle*,Particle_List> incs;
    size_t n(0), s_nmax(100000);
    for (bool cont(true);cont&&++n<s_nmax;) {
      cont=false;
      msg_Debugging()<<"  clustering "<<n<<" {\n";
      for (size_t i(0);i<outlist.size();++i) {
	if (outlist[i]->Info()=='S') continue;
	incs[outlist[i]].clear();
	for (Particle_List::iterator pit(stowlist.begin());
	     pit!=stowlist.end();++pit)
	  if (DeltaR2(outlist[i],*pit)<m_deltar2) 
	    incs[outlist[i]].push_back(*pit);
	Vec4D center(ETCenter(incs[outlist[i]]));
	msg_Debugging()<<"    new cluster "<<i<<": "
		       <<outlist[i]->Momentum()<<" -> "<<center<<" ("
		       <<sqrt(DeltaR2(outlist[i]->Momentum(),center))
		       <<" vs. "<<sqrt(m_deltarmin2)<<") {\n";
	for (Particle_List::iterator pit(incs[outlist[i]].begin());
	     pit!=incs[outlist[i]].end();++pit) 
	  msg_Debugging()<<"      "<<(*pit)->Momentum()<<"\n";
	msg_Debugging()<<"    }\n";
	if (DeltaR2(outlist[i]->Momentum(),center)<m_deltarmin2)
	  outlist[i]->SetInfo('S');
	else cont=true;
	outlist[i]->SetMomentum(center);
      }
      msg_Debugging()<<"  }\n";
    }
    if (n==s_nmax) 
      msg_Error()<<METHOD<<"(): Caught in infinite loop. Abort."<<std::endl;
    // set cluster momenta
    for (size_t i(0);i<outlist.size();++i) {
      msg_Debugging()<<"  set momentum "<<i<<": "<<outlist[i]->Momentum();
      outlist[i]->SetMomentum(ECenter(incs[outlist[i]]));
      msg_Debugging()<<" -> "<<outlist[i]->Momentum()
		     <<" "<<outlist[i]<<"\n";
    }
    // sort in et
    std::stable_sort(outlist.begin(),outlist.end(),Order_ET());
    // merge clusters
    msg_Debugging()<<"  merging jets {\n";
    for (bool cont(true);cont;) {
      cont=false;
      for (size_t i(0);i<outlist.size()-1;++i) {
	if (outlist[i]==NULL) continue;
	msg_Debugging()<<"    check "<<i<<std::flush<<": "<<outlist[i]
		       <<" { "<<outlist[i]->Momentum()<<" {\n";
	for (Particle_List::iterator pit(incs[outlist[i]].begin());
	     pit!=incs[outlist[i]].end();++pit) 
	  msg_Debugging()<<"        "<<(*pit)->Momentum()<<"\n";
	msg_Debugging()<<"      }\n";
	for (size_t j(i+1);j<outlist.size();++j) {
	  if (outlist[j]==NULL) continue;
	  msg_Debugging()<<"      check "<<j<<std::flush<<": "<<outlist[j]
			 <<" "<<outlist[j]->Momentum()<<" {\n";
	  for (Particle_List::iterator pit(incs[outlist[j]].begin());
	       pit!=incs[outlist[j]].end();++pit) 
	    msg_Debugging()<<"          "<<(*pit)->Momentum()<<"\n";
	  msg_Debugging()<<"        }\n";
	  Particle_List tlist;
	  for (Particle_List::iterator pit(incs[outlist[i]].begin());
	       pit!=incs[outlist[i]].end();++pit) {
	    if (DeltaR2(outlist[j],*pit)<m_deltar2) tlist.push_back(*pit);
	  }
	  if (tlist.size()) {
	    msg_Debugging()<<"      -> "<<2.0*ECenter(tlist).EPerp()
			   <<" vs. Min("<<outlist[i]->Momentum().EPerp()
			   <<","<<outlist[j]->Momentum().EPerp()<<")\n";
	    if (ECenter(tlist).EPerp()<0.5*Min
		(outlist[i]->Momentum().EPerp(),
		 outlist[j]->Momentum().EPerp())) {
	      msg_Debugging()<<"        reassigning "<<i<<" & "<<j;
	      bool hit(false);
	      for (Particle_List::iterator pit(tlist.begin());
		   pit!=tlist.end();++pit)
		if (DeltaR2(outlist[i],*pit)<DeltaR2(outlist[j],*pit)) {
		  Add(incs[outlist[i]],*pit);
		  if (Remove(incs[outlist[j]],*pit)) hit=true;
		}
		else {
		  Add(incs[outlist[j]],*pit);
		  if (Remove(incs[outlist[i]],*pit)) hit=true;
		}
	      if (hit) cont=true;
	      msg_Debugging()<<" -> "<<hit<<"\n";
	      outlist[i]->SetMomentum(ECenter(incs[outlist[i]]));
	      outlist[j]->SetMomentum(ECenter(incs[outlist[j]]));
	    }
	    else {
	      msg_Debugging()<<"        merging "<<i<<" & "<<j<<"\n";
	      for (Particle_List::iterator pit(tlist.begin());
		   pit!=tlist.end();++pit) Add(incs[outlist[i]],*pit);
	      outlist[i]->SetMomentum(ECenter(incs[outlist[i]]));
	      delete outlist[j];
	      outlist[j]=NULL;
 	      cont=true;
	    }
	    if (cont) {
	      std::stable_sort(outlist.begin(),outlist.end(),Order_SET());
	      i=0;
	      j=0;
	      cont=false;
	      msg_Debugging()<<"    }\n";
	      msg_Debugging()<<"    check "<<i<<std::flush<<": "<<outlist[i]
			     <<" { "<<outlist[i]->Momentum()<<" {\n";
	      for (Particle_List::iterator pit(incs[outlist[i]].begin());
		   pit!=incs[outlist[i]].end();++pit) 
		msg_Debugging()<<"        "<<(*pit)->Momentum()<<"\n";
	      msg_Debugging()<<"      }\n";
	    }
	  }
	}
	msg_Debugging()<<"    }\n";
      }
    }
    msg_Debugging()<<"  }\n";
    // add final jet list
    for (Particle_List::iterator pit(outlist.begin());
 	 pit!=outlist.end();) {
      if (*pit==NULL) outlist.erase(pit);
      else ++pit;
    }
    msg_Debugging()<<"  final jet list{\n";
    for (Particle_List::iterator pit(outlist.begin());
 	 pit!=outlist.end();++pit) {
      msg_Debugging()<<"    "<<**pit<<"\n";
      for (Particle_List::iterator spit(incs[*pit].begin());
	   spit!=incs[*pit].end();++spit) 
	msg_Debugging()<<"      "<<**spit<<"\n";
    }
    msg_Debugging()<<"  }\n";
    msg_Debugging()<<"}\n";
  }

};

DECLARE_GETTER(DB_Cone_Getter,"DBCone",Analysis_Object,Argument_Matrix);

Analysis_Object *DB_Cone_Getter::
operator()(const Argument_Matrix &parameters) const	
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<6) return NULL;
    return new 
      D0_BFKL_Cone_Finder(ToType<double>(parameters[0][0]),
			  ToType<double>(parameters[0][1]),
			  ToType<double>(parameters[0][2]),
			  ToType<double>(parameters[0][3]),
			  parameters[0][4],parameters[0][5]);
  }
  else if (parameters.size()<6) return NULL;
  double etmin=1.0, deltarpre=0.3, deltar=0.7, deltarmin=0.001;
  std::string inlist="FinalState", outlist="LeadPartJets";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="ETMin") etmin=ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="DeltaRPre") deltarpre=ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="DeltaR") deltar=ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="DeltaRMin") deltarmin=ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
  }
  return new 
    D0_BFKL_Cone_Finder(etmin,deltarpre,deltar,deltarmin,inlist,outlist);
}									

void DB_Cone_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"etmin deltarpre deltar deltarmin inlist outlist"; 
}


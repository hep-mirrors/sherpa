#include "QCD_Remnant_Base.H"

#include "Exception.H"
#include "Random.H"
#include <iomanip>

#ifdef PROFILE__all
#define PROFILE__QCD_Remnant_Base
#endif
#ifdef PROFILE__QCD_Remnant_Base
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

QCD_Remnant_Base::QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,const unsigned int beam,
				   const double scale,const rtp::code type):
  Remnant_Base(type,beam),
  m_deltax(0.0125), 
  m_scale(scale),
  m_ecms(sqrt(isrhandler->Pole())),
  m_xscheme(1), 
  m_maxtrials(100)
{
  if (isrhandler==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"QCD remnant needs ISR Handler.",
			    "QCD_Remnant_Base","QCD_Remnant_Base"));
  }
  p_pdfbase=isrhandler->PDF(m_beam)->GetBasicPDF();
  m_dupdf=isrhandler->KMROn()>0;
  m_finder.Exclude(ATOOLS::btp::Beam);
}

void QCD_Remnant_Base::Clear()
{
  while (m_sorted.size()>0) {
    delete m_sorted.front();
    m_sorted.erase(m_sorted.begin());
  }
  Remnant_Base::Clear();
}

bool QCD_Remnant_Base::TestColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
				 bool singlet,bool force,int anti)
{
  if (particle->GetFlow(anti)==oldc && 
      (m_adjusted.find(particle)==m_adjusted.end() || 
       (force && m_singlet.find(particle)==m_singlet.end()))) {
    return true;
  }
  return false;
}

bool QCD_Remnant_Base::AdjustColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
				   bool &singlet,bool force,int anti,bool forward)
{
  if (m_adjusted.find(particle)!=m_adjusted.end()) m_singlet.insert(particle);
  m_adjusted.insert(particle);
  if (!force && (particle->GetFlow(1)==newc || particle->GetFlow(2)==newc)) {
    msg_Tracking()<<"QCD_Remnant_Base::AdjustColours(..): "
		  <<"Created colour singlet. Retry."<<std::endl;
    return singlet=true;
  }
  msg_Debugging()<<*particle<<" "<<oldc<<" -> "<<newc<<" "<<anti<<" "<<forward<<std::endl;
  if ((forward && particle->DecayBlob()==NULL) ||
      (!forward && particle->ProductionBlob()==NULL)) {
    return true;
  }
  if (m_adjusted.size()>100) {
    ATOOLS::msg.Error()<<"QCD_Remnant_Base::AdjustColours(..): "
		       <<"Colour nesting is too deep (more than "<<m_adjusted.size()-1
		       <<" levels)."<<std::endl
		       <<"   Cannot adjust colours completely. "
		       <<"Result might be unreliable."<<std::endl;
    return false;
  }
  ATOOLS::Blob *cur=particle->DecayBlob();
  int newanti=anti;
  bool newforward=forward;
  if (forward) {
    for (int i=0;i<cur->NOutP();++i) {
      ATOOLS::Particle *help=cur->OutParticle(i);
      if (TestColours(help,oldc,newc,singlet,force,newanti)) { 
	if (!AdjustColours(help,oldc,newc,singlet,force,newanti,newforward)) return false;
	if (!singlet) help->SetFlow(newanti,newc);
	return true;
      }
    }
    newanti=3-newanti;
    newforward=!newforward;
    for (int i=0;i<cur->NInP();++i) {
      ATOOLS::Particle *help=cur->InParticle(i);
      if (TestColours(help,oldc,newc,singlet,force,newanti)) { 
	if (!AdjustColours(help,oldc,newc,singlet,force,newanti,newforward)) return false;
	if (!singlet) help->SetFlow(newanti,newc);
	return true;
      }
    }
  }
  else {
    cur=particle->ProductionBlob();
    for (int i=0;i<cur->NInP();++i) {
      ATOOLS::Particle *help=cur->InParticle(i);
      if (TestColours(help,oldc,newc,singlet,force,newanti)) { 
	if (!AdjustColours(help,oldc,newc,singlet,force,newanti,newforward)) return false;
	if (!singlet) help->SetFlow(newanti,newc);
	return true;
      }
    }
    newanti=3-newanti;
    newforward=!newforward;
    for (int i=0;i<cur->NOutP();++i) {
      ATOOLS::Particle *help=cur->OutParticle(i);
      if (TestColours(help,oldc,newc,singlet,force,newanti)) { 
	if (!AdjustColours(help,oldc,newc,singlet,force,newanti,newforward)) return false;
	if (!singlet) help->SetFlow(newanti,newc);
	return true;
      }
    }
  }
  return true;
}

bool QCD_Remnant_Base::AdjustColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
				   bool &singlet,bool force)
{
  PROFILE_HERE;
  m_singlet.clear();
  m_adjusted.clear();
  if (oldc==newc) return true; 
  size_t i=1;
  for (;i<3;++i) if (particle->GetFlow(i)==oldc) break;
  ATOOLS::Parton_Finder finder;
  particle=finder.FindConnected(particle,true,i);
  for (i=1;i<3;++i) if (particle->GetFlow(i)==oldc) break;
  bool result=AdjustColours(particle,oldc,newc,singlet,force,i,particle->DecayBlob()!=NULL);
  if (result && !singlet) { 
    particle->SetFlow(i,newc);
  }
  return result;
}

ATOOLS::Particle *QCD_Remnant_Base::FindConnected(ATOOLS::Particle *final,
						  const bool anti) 
{
  PROFILE_HERE;
  ATOOLS::Particle *selected=NULL;
  for (short unsigned int set=0;set<2;++set) {
    for (size_t i=0;i<m_parton[set].size();++i) {
      if (m_parton[set][i]->GetFlow(2-anti)==final->GetFlow(1+anti)) {
	if (m_adjusted.find(m_parton[set][i])==m_adjusted.end()) {
	  selected=m_parton[set][i];
	}
      }
    }
  }
  return selected;
}

ATOOLS::Particle *QCD_Remnant_Base::FindDisconnected(ATOOLS::Particle *final,
						     const bool anti) 
{
  PROFILE_HERE;
  ATOOLS::Particle *selected=NULL;
  for (short unsigned int set=0;set<2;++set) {
    for (size_t i=0;i<m_parton[set].size();++i) {
      if (m_adjusted.find(m_parton[set][i])==m_adjusted.end()) {
	msg_Debugging()<<"take "<<m_parton[set][i]<<" "<<" "<<final<<std::endl;
	selected=m_parton[set][i];
      }
      else if (m_initial==2) {
	if (m_finder.FindConnected(m_parton[set][i])!=final) selected=m_parton[set][i];
      }
    }
  }
  return selected;
}

ATOOLS::Particle *QCD_Remnant_Base::FindClosest(ATOOLS::Particle *final,
						const bool anti) 
{
  PROFILE_HERE;
  ATOOLS::Particle *selected=NULL;
  double min=std::numeric_limits<double>::max();
  for (short unsigned int set=0;set<2;++set) {
    for (size_t i=0;i<m_parton[set].size();++i) {
      if (m_adjusted.find(m_parton[set][i])==m_adjusted.end()) {
	ATOOLS::Particle *end=m_finder.FindConnected(m_parton[set][i]);
	double dist=end->Momentum().PPerp(final->Momentum());
	if (min>dist && dist>0. && m_parton[set][i]->GetFlow(2-anti)!=0) {
	  msg_Debugging()<<"compare "<<m_parton[set][i]<<" "<<end<<" "<<" "<<final<<std::endl;
	  min=dist;
	  selected=m_parton[set][i];
	}
      }
    }
  }
  return selected;
}

bool QCD_Remnant_Base::SelectCompanion(QCD_Remnant_Info *const cur) 
{
  PROFILE_HERE;
  msg_Debugging()<<"look for "<<*cur<<std::endl;
  m_finder.SetColour(1+!*cur,(*cur)->GetFlow(1+!cur));
  ++*cur=m_finder.FindConnected((*cur)(),true);
  m_adjusted.clear();
  size_t k=0;//, l=0;
  for (short unsigned int i=0;i<2;++i) {
    m_adjusted.insert((*cur)(1));
    while ((*cur)[i]==NULL) {
      //      msg_Debugging()<<++l<<std::endl;
      ATOOLS::Particle *final=FindConnected((*cur)(i),i);
      if (final==NULL) {
	final=FindClosest(++*cur,i);
	if (final==NULL) {
	  final=FindDisconnected(++*cur,i);
	  if (final==NULL) {
	    (*cur)[i]=cur;
	    continue;
	  }
	}
      }
      m_finder.SetColour(2-i,final->GetFlow(2-i));
      ATOOLS::Particle *initial=m_finder.FindConnected(final,false);
      m_adjusted.insert(initial);
      m_adjusted.insert(final);
      //      msg_Debugging()<<*initial<<" "<<initial<<" "<<final<<" "
      //			     <<++*cur<<" "<<m_adjusted.size()<<" ("<<m_constrained<<" "
      //			     <<p_beamblob->NOutP()<<")"<<std::endl;
      if (m_constrained<m_initial-1) if (cur->Find(initial)>-1) continue;
      //      msg_Debugging()<<"constrained "<<m_constrained<<" vs "<<m_initial<<" "<<(cur->Find(initial)>-1)<<std::endl;
      //      msg_Debugging()<<*initial<<" "<<initial<<" "<<final<<" "
      //                             <<++*cur<<" "<<m_adjusted.size()<<std::endl;
//       for (std::set<ATOOLS::Particle*>::iterator pit=m_adjusted.begin();pit!=m_adjusted.end();++pit) 
// 	msg_Debugging()<<*pit<<" "; msg_Debugging()<<(m_adjusted.find(initial)==m_adjusted.end())<<std::endl;
      for (size_t j=0;j<m_sorted.size();++j) {
	if ((*m_sorted[j])(0)==initial || (*m_sorted[j])(1)==initial) {
	  ++k;
//  	  msg_Debugging()<<i<<" "<<++k<<" "<<m_adjusted.size()<<" "<<m_sorted.size()<<" "
// 				 <<m_parton[1].size()<<" found "<<*m_sorted[j]<<" "
// 				 <<((*m_sorted[j])[1-i]!=NULL)<<std::endl;
	  if ((*m_sorted[j])[1-i]!=NULL) {
	    break;
	  }
	  else {
	    (*cur)[i]=m_sorted[j];
	    (*m_sorted[j])[1-i]=cur;
	    m_adjusted.clear();
	    m_adjusted.insert((*cur)(0));
	    m_adjusted.insert(final);
	    ++m_constrained;
	    break;
	  }
	}
      }
    }
  }
  //  msg_Debugging()<<"leave "<<*cur<<"__________________________________________________"<<std::endl;
  return true;
}

void QCD_Remnant_Base::SortRemnants() 
{
  PROFILE_HERE;
  m_sorted.clear();
  for (size_t i=0;i<m_parton[1].size();++i) {
    QCD_Remnant_Info *shifted = new QCD_Remnant_Info(m_parton[1][i]);
    if ((*shifted)(0)!=(*shifted)(1)) {
      ATOOLS::Particle *comp=(*shifted)(1);
      if ((*shifted)(1)==(*shifted)()) comp=(*shifted)(0);
      p_beamblob->AddToOutParticles(comp);
      m_parton[0].push_back(comp);
    }
    for (size_t j=0;j<m_sorted.size();++j) {
      if (m_finder.FindConnected((*m_sorted[j])())->Momentum().PPerp2()<
	  m_finder.FindConnected(m_parton[1][i])->Momentum().PPerp2()) {
	std::swap<QCD_Remnant_Info*>(shifted,m_sorted[j]);
	for (++j;j<m_sorted.size();++j) {
	  std::swap<QCD_Remnant_Info*>(shifted,m_sorted[j]);
	}
      }
    }
    m_sorted.push_back(shifted);
  }
  m_initial=m_sorted.size();
//   for (size_t j=0;j<m_sorted.size();++j) {
//     msg_Debugging()<<j<<" -> "<<m_finder.FindConnected((*m_sorted)[j]())
//       ->Momentum().PPerp2()<<*(*m_sorted)[j]()<<"\n";
//   }
}

bool QCD_Remnant_Base::SelectCompanions() 
{
  PROFILE_HERE;
  m_constrained=0;
  for (size_t i=0;i<m_sorted.size();++i) {
    if (!SelectCompanion(m_sorted[i])) return false; 
  }
  return true;
}

bool QCD_Remnant_Base::ConnectRemnants()
{
  PROFILE_HERE;
  m_adjusted.clear();
  for (size_t i=0;i<m_sorted.size();++i) {
    for (size_t j=0;j<2;++j) {
      if (!TreatDipole(m_sorted[i],j)) {
	UnDo();
	p_partner->UnDo();
	return false;
      }
    }
  }
//   msg_Debugging()<<"=====next===="<<std::endl;
//   msg_Debugging()<<*p_beamblob<<std::endl;
  return true;
}

bool QCD_Remnant_Base::TreatDipole(QCD_Remnant_Info *const cur,const size_t i) 
{
  PROFILE_HERE;
  //  msg_Debugging()<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"<<*cur<<std::endl;
  if (m_adjusted.find((*cur)(i))!=m_adjusted.end()) return true;
  m_adjusted.insert((*(*cur)[i])(1-i));
  bool singlet=false;
  unsigned int old=(*(*cur)[i])(1-i)->GetFlow(2-i);
  if (!AdjustColours((*(*cur)[i])(1-i),old,(*cur)(i)->GetFlow(1+i),singlet,false)) return false;
  if ((*(*cur)[i])(i)->GetFlow(1+i)==(*cur)->GetFlow(2-i)) {
    p_beamblob->RemoveOutParticle((*(*cur)[i])(1-i));
    AdjustColours((*(*cur)[i])(1-i),(*cur)(i)->GetFlow(1+i),old,singlet,false);
    p_beamblob->AddToOutParticles((*(*cur)[i])(1-i));
    singlet=true;
  }
  if (singlet && (*(*cur)[i])(1-i)->GetFlow(2-i)!=(*cur)(i)->GetFlow(1+i)) {
    ATOOLS::Particle *newpart = new ATOOLS::Particle(-1,ATOOLS::kf::gluon);
    QCD_Remnant_Info *newinfo = new QCD_Remnant_Info(newpart);
    m_parton[0].push_back(newpart);
    m_sorted.push_back(newinfo);
    newpart->SetFlow(2-i,(*cur)(i)->GetFlow(1+i));
    newpart->SetFlow(1+i,(*(*cur)[i])(1-i)->GetFlow(2-i));
    p_beamblob->AddToOutParticles(newpart);
    (*newinfo)[1-i]=cur;
    (*newinfo)[i]=(*cur)[i];
    (*(*cur)[i])[1-i]=newinfo;
    (*cur)[i]=newinfo;
    //    msg_Debugging()<<"inserted "<<i<<" "<<*newinfo<<std::endl;
  }
  else (*cur)-i=old;
  //  msg_Debugging()<<*p_beamblob<<std::endl;
  return true;
}

void QCD_Remnant_Base::UnDo() 
{
  PROFILE_HERE;
  msg_Tracking()<<"QCD_Remnant_Base::UnDo(): Undoing changes on blob list."<<std::endl;
  for (int i=(int)m_sorted.size()-1;i>=0;--i) {
    for (int j=1;j>=0;--j) {
      bool singlet=false;
      AdjustColours((*(*m_sorted[i])[j])(1-j),
 		    (*(*m_sorted[i])[j])(1-j)->GetFlow(2-j),(*m_sorted[i])-j,singlet,true);
    }
  }
  while (m_sorted.size()>0) {
    delete m_sorted.front();
    m_sorted.erase(m_sorted.begin());
  }
  Remnant_Base::UnDo();
}


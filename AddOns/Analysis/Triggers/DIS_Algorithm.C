#include "AddOns/Analysis/Triggers/DIS_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

#include <iomanip>
#include <algorithm>
#include <limits>

#ifdef PROFILE__all
#define PROFILE__Analysis_Phase
#endif
#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

DIS_Algorithm::DIS_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier) : 
  Jet_Algorithm_Base(qualifier), m_matrixsize(0), p_ktij(NULL), p_imap(NULL), 
  p_kis(NULL), p_jets(NULL), p_kts(NULL)
{ 
}

DIS_Algorithm::~DIS_Algorithm()
{
  if (p_ktij) {
    for (int i=0;i<m_matrixsize;++i) {
      delete [] p_ktij[i];
    }
    delete [] p_ktij;
    delete [] p_imap;
    delete [] p_kis;
  }
}


void DIS_Algorithm::AddToKtlist(double kt2) {
  if (p_kts) {
    p_kts->push_back(kt2);
  }
}

void DIS_Algorithm::AddToJetlist(const Vec4D & mom, int bf) {
  if (p_jets) {
    if(!bf) p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_jet),mom));
    else p_jets->push_back(new Particle(p_jets->size(),bf>0?Flavour(kf_bjet):
					Flavour(kf_bjet).Bar(),mom));
  }
}



bool DIS_Algorithm::ConstructJets(const Particle_List * pl, Particle_List * jets, 
				  std::vector<double> * kts, double rmin) {
  // assume empty containers
  p_jets = jets;
  p_kts  = kts;
  m_r2min = sqr(rmin);

  // create vector list from particle list
  
  int n=0;
  Vec4D * moms = new Vec4D[pl->size()];
  int * bflag = new int[pl->size()];
  static Particle_Qualifier_Base 
    *bhad(Particle_Qualifier_Getter::GetObject("DecayedBHadron",""));
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    if ((*p_qualifier)(*it)) {
      //    if (!(*it)->Flav().IsLepton()) {
      moms[n]  = ((*it)->Momentum()); 
      bflag[n] = false;
      if (m_bflag==0) bflag[n] = (((*it)->Flav()).Kfcode()==kf_b || 
				  (*bhad)(*it) ||
				  ((*it)->Flav()).Kfcode()==kf_bjet);
      else if (m_bflag==-1) bflag[n] = (1-2*(*it)->Flav().IsAnti())*
			      (((*it)->Flav()).Kfcode()==kf_b || 
			       (*bhad)(*it) ||
			       ((*it)->Flav()).Kfcode()==kf_bjet);
      ++n;
    }
  }

  // boost into breit frame
  Vec4D l, lp, pp;
  for (size_t i(0);i<p_bl->size();++i)
    if ((*p_bl)[i]->Type()==btp::Beam) {
      Particle *p((*p_bl)[i]->InParticle(0));
      if (p->Flav().IsLepton()) l=p->Momentum();
      else pp=p->Momentum();
    }
  Blob *me(p_bl->FindFirst(btp::Signal_Process));
  for (int i(0);i<me->NOutP();++i)
    if (me->OutParticle(i)->Flav().IsLepton()) {
      lp=me->OutParticle(i)->Momentum();
      break;
    }
  Vec4D qq(l-lp);
  Poincare cms(pp+qq);
  double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), E(sqrt(Q2)/(2.0*x));
  m_p=Vec4D(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
  m_q=Vec4D(0.0,0.0,0.0,2.0*x*E);
  cms.Boost(pp);
  cms.Boost(qq);
  Poincare zrot(pp,-Vec4D::ZVEC);
  zrot.Rotate(pp);
  zrot.Rotate(qq);
  Poincare breit(m_p+m_q);
  breit.BoostBack(pp);
  breit.BoostBack(qq);
  if (!IsEqual(pp,m_p,1.0e-3) || !IsEqual(qq,m_q,1.0e-3))
    msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
  for (int i(0);i<n;++i) {
    cms.Boost(moms[i]);
    zrot.Rotate(moms[i]);
    breit.BoostBack(moms[i]);
  }

  // cluster
  Ktmin(moms,bflag,n);
  delete [] moms;
  delete [] bflag;

  // finalize (sort and release used containers)
  SortPT(p_jets);

  for (size_t i(0);i<p_jets->size();++i) {
    Vec4D p((*p_jets)[i]->Momentum());
    breit.Boost(p);
    zrot.RotateBack(p);
    cms.BoostBack(p);
    (*p_jets)[i]->SetMomentum(p);
  }

  p_jets=0;
  p_kts =0;

  return true;
}


void DIS_Algorithm::Init(int size) 
{
  PROFILE_HERE;
  if (size>m_matrixsize ) {
    if (p_ktij) {
      for (int i=0;i<m_matrixsize;++i) delete [] p_ktij[i];
      delete [] p_ktij;
      delete [] p_imap;
      delete [] p_kis;
    }
    m_matrixsize = size;
    p_kis  = new double[size];
    p_imap = new int[size];
    p_ktij = new double*[size];
    for (int i=0;i<size;++i) p_ktij[i] = new double[size];
  }
  for (int i=0;i<size;++i) p_imap[i] = i;
}


double DIS_Algorithm::Ktmin(Vec4D * p, int * bf, int n)
{
  if (n==0) return 0.;
  if (n==1) {
    AddToJetlist(p[0],bf[0]);
    double dmin=sqr(p[0][0])*R2(p[0]);
    AddToKtlist(dmin);
    return dmin;;
  }

  Init(n);


  //cal first matrix
  int ii=0, jj=0;
  double dmin=std::numeric_limits<double>::max();
  {
    PROFILE_LOCAL(" first loop ");
    
    for (int i=0;i<n;++i) {
      double di = p_ktij[i][i] = sqr(p[i][0])*R2(p[i]);
      if (di<dmin) { dmin=di; ii=i; jj=i;}
      for (int j=0;j<i;++j) {
	double dij = p_ktij[i][j] = sqr(Min(p[i][0],p[j][0]))*R2(p[i],p[j]);
	if (dij<dmin) {dmin=dij; ii=i; jj=j;}
      }
    }
  }

  // recalc matrix
  while (n>0) {
    PROFILE_LOCAL(" main loop ");
    if (dmin>m_r2min*sqr(rpa.gen.Ecms())) break;
    if (ii!=jj) {
      // combine precluster
#ifdef DEBUG_JETRATES
      msg_Debugging()<<"jetrate Q_{"<<n<<"->"
		     <<n-1<<"} = "<<sqrt(dmin)
		     <<" <- "<<p[p_imap[jj]]<<" "<<p[p_imap[ii]]<<"\n";
#endif
      p[p_imap[jj]]+=p[p_imap[ii]];
      bf[p_imap[jj]] = bf[p_imap[jj]]+bf[p_imap[ii]];      
      AddToKtlist(dmin);
    }
    else {
      // add to jet list
#ifdef DEBUG_JETRATES
      msg_Debugging()<<"jetrate Q_{"<<n<<"->"
		     <<n-1<<"} = "<<sqrt(dmin)
		     <<" <- "<<p[p_imap[jj]]<<"\n";
#endif
      AddToKtlist(dmin);
    }

    --n;

    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n==1) {
      int jjx=p_imap[jj];
      p_ktij[jjx][jjx] = sqr(p[jjx][0])*R2(p[jjx]);
      break;
    }
    // update map (remove precluster)
    {
      PROFILE_LOCAL(" update loop ");

    
    // update matrix (only what is necessary)
    int jjx=p_imap[jj];
    p_ktij[jjx][jjx] = sqr(p[jjx][0])*R2(p[jjx]);
    for (int j=0;j<jj;++j)   p_ktij[jjx][p_imap[j]] = 
      sqr(Min(p[jjx][0],p[p_imap[j]][0]))
      *R2(p[jjx],p[p_imap[j]]);
    for (int i=jj+1;i<n;++i) p_ktij[p_imap[i]][jjx] = 
      sqr(Min(p[jjx][0],p[p_imap[i]][0]))
      *R2(p[jjx],p[p_imap[i]]);
    }
    // redetermine rmin and dmin
    ii=jj=0;
    {
      PROFILE_LOCAL(" second loop ");

    dmin=p_ktij[p_imap[0]][p_imap[0]];
    for (int i=0;i<n;++i) {
      int ix=p_imap[i];
      double di = p_ktij[ix][ix];
      if (di<dmin) { dmin=di; ii=jj=i;}
      for (int j=0;j<i;++j) {
	int jx=p_imap[j];
	double dij = p_ktij[ix][jx];
	if (dij<dmin) { dmin=dij; ii=i; jj=j;}
      }
    }
    }
  }

  // add remaining preclusters to jetlist
  for (int i=0;i<n;++i) {
    AddToJetlist(p[p_imap[i]],bf[p_imap[i]]);
#ifdef DEBUG_JETRATES
    msg_Debugging()<<"jet at rate Q_{"<<n<<"->"
		   <<n-1<<"} = "<<sqrt(p_ktij[p_imap[i]][p_imap[i]])
		   <<" <- "<<p[p_imap[i]]<<"\n";
#endif
    AddToKtlist(p_ktij[p_imap[i]][p_imap[i]]);
  }
  return dmin;
}

double DIS_Algorithm::R2(const Vec4D &p1) const
{
  return 2.0*(1.0-p1.CosTheta(m_p));
}

double DIS_Algorithm::R2(const Vec4D &p1, const Vec4D &p2) const
{
  return 2.0*(1.0-p1.CosTheta(p2));
}


#include "AddOns/Analysis/Triggers/Cambridge_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Message.H"

#include <iomanip>
#include <algorithm>

//#define PROFILE__Analysis_Phase

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;


Cambridge_Algorithm::Cambridge_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier, int mode) : 
 Jet_Algorithm_Base(qualifier),  m_mode(mode), m_matrixsize(0), p_yij(NULL), p_imap(NULL), m_vectorsize(0), p_moms(NULL), 
  p_bflag(NULL), p_jets(NULL)
{

}

Cambridge_Algorithm::~Cambridge_Algorithm()
{
  if (p_moms)  delete [] p_moms;
  if (p_bflag) delete [] p_bflag;

  if (p_imap) {
    delete [] p_imap;
  }
  if (p_yij) {
    for (int i=0;i<m_matrixsize;++i) delete [] p_yij[i];
    delete [] p_yij;
  }
}


void Cambridge_Algorithm::AddToKtlist(double kt2) {
   if (p_kts) {
     p_kts->push_back(sqr(kt2));
   }
}

void Cambridge_Algorithm::AddToJetlist(const Vec4D & mom, bool bf) {
  ++m_njets;
  if (p_jets) {
    if(!bf) p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_jet),mom));
    else    p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_bjet),mom));
  }
}



bool Cambridge_Algorithm::ConstructJets(const Particle_List * pl, Particle_List * jets, 
				     std::vector<double> * kts, double ycut) {
  PROFILE_HERE;
  p_jets = jets;
  p_kts  = kts;
  m_ycut = ycut;

  Vec4D momsum=Vec4D(0.,0.,0.,0.);
  // create vector list from particle list
  int n=0;

  InitMoms(pl->size());
    
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    momsum+=(*it)->Momentum();
    if ((*p_qualifier)(*it)) {
      p_moms[n]  = ((*it)->Momentum()); 
      p_bflag[n] = (((*it)->Flav()).Kfcode()==kf_b)&& !m_bflag;
      ++n;
    }
  }
  m_sprime = momsum.Abs2();
  
  if (m_mode==0) {    // generate jetlist only!
    p_kts =0;

    // cluster
    Ymin(p_moms,p_bflag,n);
    // finalize (sort and release used containers)

    SortE(p_jets);

    p_jets=0;
  
    return true;
  }
  // (m_mode==1) generate yij lists !


  m_ymax=1.;
  p_kts->push_back(1.);
  m_jetnumbers.clear();
  for (;;)  {
    p_kts =0;
    m_ycut=m_ymax-1.e-10;
    m_njets=0;
    
    m_ymax=0.;
    Ymin(p_moms,p_bflag,n);

    // store ymax and number of jets
    kts->push_back(sqr(m_ymax));
    m_jetnumbers.push_back(m_njets);

    //    std::cout<<" ymax="<<m_ymax<<"  njets="<<m_njets<<"\n";
    if (m_ymax<ycut) break;
    n=0;
    for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
      momsum+=(*it)->Momentum();
      if ((*p_qualifier)(*it)) {
	p_moms[n]  = ((*it)->Momentum()); 
	p_bflag[n] = (((*it)->Flav()).Kfcode()==kf_b)&& !m_bflag;
	++n;
      }
    }
  }
  // store information
  std::sort(kts->begin(), kts->end());

  // ...
  return true;
}


void Cambridge_Algorithm::InitMoms(int size) 
{
  if (2*size>m_vectorsize ) {
    m_vectorsize=2*size;
    if (p_moms)  delete [] p_moms;
    if (p_bflag) delete [] p_bflag;
    p_moms = new Vec4D[size*2];
    p_bflag = new bool[size*2];
  }
}

void Cambridge_Algorithm::Init(int size) 
{
  PROFILE_HERE;
  if (size>m_matrixsize ) {
    if (p_yij) {
      for (int i=0;i<m_matrixsize;++i) delete [] p_yij[i];
      delete [] p_yij;
      delete [] p_imap;
    }
    m_matrixsize = size;
    p_imap = new int[size];
    p_yij = new double*[size];
    for (int i=0;i<size;++i) p_yij[i] = new double[size];
  }
  for (int i=0;i<size;++i) p_imap[i] = i;
}


void Cambridge_Algorithm::Ymin(Vec4D * p, bool * bf, int n)
{
  PROFILE_HERE;
  if (n==0) return;
  if (n==1) {
    AddToJetlist(p[0],bf[0]);
    return;
  }

  Init(n);

  //cal first matrix
  int ii=0, jj=0;
//   double yij=0;
//   double ymin=1.;
  double dmin=4.;
  {
    //    PROFILE_LOCAL(" first loop ");
    
    for (int i=1;i<n;++i) {
      for (int j=0;j<i;++j) {
	double d = p_yij[i][j] = R12(p[i],p[j]);
	//	double y = sqr(Min(p[i][0],p[j][0]))*d;
	if (d<dmin) { dmin=d; ii=i; jj=j;}
      }
    }
  }

  // recalc matrix
//   int hit = 0;
  while (n>1) {
    //    PROFILE_LOCAL(" main loop ");
    
    double yij = sqr(Min(p[p_imap[ii]][0],p[p_imap[jj]][0]))*dmin;
    //    std::cout<<"winner : "<<ii<<","<<jj<<"  d="<<dmin<<"  y="<<yij<<"\n";
//     if (!hit && ymin>=m_ycut) {
//       hit = 1;
//       for (int i=0;i<n;++i) AddToJetlist(p[p_imap[i]],bf[p_imap[i]]);
//     }

    AddToKtlist(yij);
    if (yij<m_ycut) {
      //      std::cout<<" combining ...\n";
      if (yij>m_ymax) {
	m_ymax=yij;
	//	std::cout<<" ymax="<<m_ymax<<"\n";
      }
      // combine precluster
      p[p_imap[jj]]+=p[p_imap[ii]];
      bf[p_imap[jj]] = bf[p_imap[jj]]||bf[p_imap[ii]];      
    }
    else {
      // soft freezing
      if (p[p_imap[ii]][0]>p[p_imap[jj]][0]) {
	int kk=ii;
	ii=jj;
	jj=kk;
	//	std::cout<<" freezing ..."<<ii<<"\n";
      }
      AddToJetlist(p[p_imap[ii]],bf[p_imap[ii]]);
    }
    --n;
    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n<=2) {
      //      std::cout<<" n="<<n<<" left\n";
      break;
    }
    // update map (remove precluster)
    {
      //      PROFILE_LOCAL(" update loop ");

      // update matrix (only what is necessary)
      int jjx=p_imap[jj];
      for (int j=0;j<jj;++j)   p_yij[jjx][p_imap[j]] = R12(p[jjx],p[p_imap[j]]);
      for (int i=jj+1;i<n;++i) p_yij[p_imap[i]][jjx] = R12(p[p_imap[i]],p[jjx]);
    }
    // redetermine dmin
    ii=jj=0;
    dmin=4.;
    {
      //      PROFILE_LOCAL(" second loop ");

      for (int i=1;i<n;++i) {
	int ix=p_imap[i];
	for (int j=0;j<i;++j) {
	  int jx=p_imap[j];
	  double d = p_yij[ix][jx];
// 	  double y = sqr(Min(p[ix][0],p[jx][0]))*d;
	  if (d<dmin) {  dmin=d; ii=i; jj=j;}
// 	  if (y<ymin)   ymin=y;
	}
      }
    }
  }
  for (int i=0;i<n;++i) AddToJetlist(p[p_imap[i]],bf[p_imap[i]]);
}

class Order_E {
public:
  int operator()(const Particle * a, const Particle * b) {
    if (a->Momentum()[0] > b->Momentum()[0]) return 1;
    return 0;
  }
};

class Order_PT {
public:
  int operator()(const Particle * a, const Particle * b) {
    if (Cambridge_Algorithm::Kt2(a->Momentum()) > Cambridge_Algorithm::Kt2(b->Momentum())) return 1;
    return 0;
  }
};



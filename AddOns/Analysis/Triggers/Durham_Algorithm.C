#include "AddOns/Analysis/Triggers/Durham_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Default_Reader.H"
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


Durham_Algorithm::Durham_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier) : 
  Jet_Algorithm_Base(qualifier), m_matrixsize(0), p_yij(NULL), p_imap(NULL), m_vectorsize(0), p_moms(NULL), 
  p_bflag(NULL), p_jets(NULL)
{

}

Durham_Algorithm::~Durham_Algorithm()
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


void Durham_Algorithm::AddToKtlist(double kt2) {
   if (p_kts) {
     p_kts->push_back(sqr(kt2));
   }
}

void Durham_Algorithm::AddToJetlist(const Vec4D & mom, int bf) {
  if (p_jets) {
    if(!bf) p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_jet),mom));
    else p_jets->push_back(new Particle(p_jets->size(),bf>0?Flavour(kf_bjet):
					Flavour(kf_bjet).Bar(),mom));
  }
}



bool Durham_Algorithm::ConstructJets(const Particle_List * pl, Particle_List * jets, 
				     std::vector<double> * kts, double ycut) {
  PROFILE_HERE;
  p_jets = jets;
  p_kts  = kts;
  m_ycut = ycut;
  

  Vec4D momsum=Vec4D(0.,0.,0.,0.);
  // create vector list from particle list
  int n=0;

  // avoid new !!!!!!!!!!!!!!!!!!!!!!!!1
  InitMoms(pl->size());

  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    momsum+=(*it)->Momentum();
    if ((*p_qualifier)(*it)) {
      p_moms[n]  = ((*it)->Momentum());
      p_bflag[n] = 0; 
      if (m_bflag==0) p_bflag[n] = (((*it)->Flav()).Kfcode()==kf_b || 
				    ((*it)->Flav()).Kfcode()==kf_bjet);
      else if (m_bflag==-1) p_bflag[n] = (1-2*(*it)->Flav().IsAnti())*
			      (((*it)->Flav()).Kfcode()==kf_b || 
			       ((*it)->Flav()).Kfcode()==kf_bjet);
      ++n;
    }
  }
  m_sprime = momsum.Abs2();

  // cluster
  Ymin(p_moms,p_bflag,n);
  // finalize (sort and release used containers)

  SortE(p_jets);

  p_jets=0;
  p_kts =0;

  return true;
}


void Durham_Algorithm::InitMoms(int size) 
{
  if (2*size>m_vectorsize ) {
    m_vectorsize=2*size;
    if (p_moms)  delete [] p_moms;
    if (p_bflag) delete [] p_bflag;
    p_moms = new Vec4D[size*2];
    p_bflag = new int[size*2];
  }
}

void Durham_Algorithm::Init(int size) 
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


void Durham_Algorithm::Ymin(Vec4D * p, int * bf, int n)
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
  double ymin=1.;
  {
    //    PROFILE_LOCAL(" first loop ");
    
    for (int i=1;i<n;++i) {
      for (int j=0;j<i;++j) {
	double y = p_yij[i][j] =Y12(p[i],p[j]);
	if (y<ymin) { ymin=y; ii=i; jj=j;}
      }
    }
  }

  // recalc matrix
  int hit = 0;
  while (n>1) {
    //    PROFILE_LOCAL(" main loop ");

    if (!hit && ymin>=m_ycut) {
      hit = 1;
      for (int i=0;i<n;++i) AddToJetlist(p[p_imap[i]],bf[p_imap[i]]);
    }

    // combine precluster
    p[p_imap[jj]]+=p[p_imap[ii]];
    bf[p_imap[jj]] = bf[p_imap[jj]]+bf[p_imap[ii]];      
    AddToKtlist(ymin);
    --n;
    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n==1) break;
    // update map (remove precluster)
    {
      //      PROFILE_LOCAL(" update loop ");

      // update matrix (only what is necessary)
      int jjx=p_imap[jj];
      for (int j=0;j<jj;++j)   p_yij[jjx][p_imap[j]] = Y12(p[jjx],p[p_imap[j]]);
      for (int i=jj+1;i<n;++i) p_yij[p_imap[i]][jjx] = Y12(p[p_imap[i]],p[jjx]);
    }
    // redetermine rmin and dmin
    ii=jj=0;
    ymin=1.;
    {
      //      PROFILE_LOCAL(" second loop ");

      for (int i=1;i<n;++i) {
	int ix=p_imap[i];
	for (int j=0;j<i;++j) {
	  int jx=p_imap[j];
	  double y = p_yij[ix][jx];
	  if (y<ymin) { ymin=y; ii=i; jj=j;}
	}
      }
    }
  }
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
    if (Durham_Algorithm::Kt2(a->Momentum()) > Durham_Algorithm::Kt2(b->Momentum())) return 1;
    return 0;
  }
};

double Durham_Algorithm::DCos12(const Vec4D & p1,const Vec4D & p2) const
{
  double s  = p1[1]*p2[1]+p1[2]*p2[2]+p1[3]*p2[3];
  double b1 = p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
  double b2 = p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
  return s/sqrt(b1*b2);
  //  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double Durham_Algorithm::Y12(const Vec4D & p1, const Vec4D & p2) const
{
  static int das(-1);
  if (das<0) {
    Default_Reader reader;
    if (reader.Read(das, "DURHAM_SCHEME", 0)) {
      msg_Info()<<METHOD<<"(): Setting durham scheme mode "<<das<<".\n";
    }
  }
  if (das==0) return 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2))/m_sprime;
  return 2.*Min(p1.PSpat2(),p2.PSpat2())*(1.-DCos12(p1,p2))/m_sprime;
}


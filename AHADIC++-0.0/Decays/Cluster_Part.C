#include "Cluster_Part.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part(const double tension,const int mode) : 
  m_popper(Pair_Popper(tension)), m_4Qmode(mode), m_ana(false) 
{ }

bool Cluster_Part::Veto(Cluster * cluster,const ATOOLS::Flavour & flav)
{
  switch (m_4Qmode) {
  case 1:
  default:
    if (flav.IsDiQuark() && 
	(cluster->GetFlav(1).IsDiQuark() || cluster->GetFlav(2).IsDiQuark()))
      return true;
  }
  return false;
}


Return_Value::code Cluster_Part::TestDecay(Cluster * const cluster)
{
  //cout<<"------------------------------------------------------------------"<<endl
  //   <<"In "<<METHOD<<endl<<(*cluster)<<endl;
  Vec4D  * momenta = new Vec4D[4];
  double * masses  = new double[4];
  double M   = cluster->Mass(0);
  masses[0]  = cluster->Mass(1);
  masses[3]  = cluster->Mass(2);
  double pz2 = (sqr(M*M-masses[0]*masses[0]-masses[3]*masses[3])-
		4.*masses[0]*masses[3])/(4.*M*M);
  if (pz2<0.) return Return_Value::Error;

  cluster->BoostInCMSAndRotateOnZ();

  double pz = sqrt(pz2); 
  double E0 = (M*M+masses[0]*masses[0]-masses[3]*masses[3])/(2.*M);
  double E3 = (M*M-masses[0]*masses[0]+masses[3]*masses[3])/(2.*M);
  momenta[0] = Vec4D(E0,0.,0.,pz); 
  momenta[3] = Vec4D(E3,0.,0.,-pz);

  Flavour flav;

  if (!BuildKinematics(cluster,flav,momenta,masses)) {
    cluster->RotateAndBoostBack();
    return Return_Value::Nothing;
  }

  Cluster * left  = new Cluster(cluster->GetFlav(1),momenta[0],flav.Bar(),momenta[1]);
  Cluster * right = new Cluster(flav,momenta[2],cluster->GetFlav(2),momenta[3]);
  if (cluster->GetLeads()==1 || cluster->GetLeads()==3) left->SetLeads(ltp::leadingtrip);
  if (cluster->GetLeads()==2 || cluster->GetLeads()==3) right->SetLeads(ltp::leadinganti);

  delete momenta;
  delete masses;

  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  cout<<METHOD<<" before boost and rotate back "<<endl<<(*cluster)<<endl;
  cluster->RotateAndBoostBack();
  cout<<METHOD<<" after boost and rotate back, boost it again "<<endl<<(*cluster)<<endl;

  Particle * part;
  part = new Particle(0,Flavour(kf::cluster),cluster->GetLeft()->Momentum()); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('C');
  part->SetFinalMass();
  cluster->GetLeft()->SetSelf(part);
  
  part = new Particle(0,Flavour(kf::cluster),cluster->GetRight()->Momentum()); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('C');
  part->SetFinalMass();
  cluster->GetRight()->SetSelf(part);

  //cout<<"Out "<<METHOD<<endl<<(*cluster)
  //   <<"------------------------------------------------------------------"<<endl;
      
  return Return_Value::Success;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Schwinger Types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Schwinger_Uniform::Schwinger_Uniform() :
  m_tension(hadpars.Get(string("Tension"))), 
  m_fraction(sqr(hadpars.Get(string("MassFraction")))),
  m_ystar_expvalue(sqr(hadpars.Get(string("<Y>")))),
  m_ystar_sigma(sqr(hadpars.Get(string("Y*_WIDTH")))),
  Cluster_Part(hadpars.Get(string("Tension")),0)
{ }

Schwinger_Uniform::~Schwinger_Uniform() { 
  cout<<METHOD<<" with analysis = "<<m_ana<<endl;
  if (m_ana) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = string("Analysis/")+hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}

void Schwinger_Uniform::SetDecayAnalysisOn() {
  m_ana = true;
  m_histograms[string("YStar_by_YStarMax")] = new Histogram(0,-1.,1.,100);
  m_histograms[string("YStar")]             = new Histogram(0,-5.,5.,100);
  m_histograms[string("YBar_by_YBarMax")]   = new Histogram(0,-1.,1.,100);
  m_histograms[string("YBar")]              = new Histogram(0,-5.,5.,100);
  m_histograms[string("Flavour")]           = new Histogram(0,0.,5.,5);
  m_histograms[string("MT2")]               = new Histogram(0,0.,5.,100);
  m_histograms[string("MT2_{u,d}")]         = new Histogram(0,0.,5.,100);
  m_histograms[string("MT2_s")]             = new Histogram(0,0.,5.,100);
  m_histograms[string("MT2_ud")]            = new Histogram(0,0.,5.,100);
  m_histograms[string("MT2_us")]            = new Histogram(0,0.,5.,100);
  m_histograms[string("MT2_ss")]            = new Histogram(0,0.,5.,100);
  m_histograms[string("SQQ")]               = new Histogram(0,0.,5.,100);
}

bool Schwinger_Uniform::BuildKinematics(Cluster * const cluster,Flavour & flav,
					Vec4D * momenta,double * masses)
{
  flav         = Flavour(kf::none);
  double MC    = cluster->Mass(0), MC2 = sqr(MC), mt2;
  double mmax2 = m_fraction*(MC2-sqr(masses[0]+masses[3]));

  double ystar, s_qq, ybmax, ybar, x1, x2, E_r, E0_r(0.), E3_r(0.);
  int counter(0);
  do {
    if (counter++>10) return false;
    switch (int(m_popper.Pop(flav,mt2,mmax2))) {
    case Return_Value::Success:
    case Return_Value::Warning:
      break;
    case Return_Value::Error: return false;
    }

    ystar      = GetYStar(mt2,mmax2);
    s_qq       = 4.*mt2*sqr(cosh(ystar));
    ybmax      = 0.5 * log(mmax2/s_qq);
    ybar       = (-1.+2.*ran.Get()) * ybmax;
    x1         = sqrt(s_qq/MC2)*exp(ybar);
    x2         = sqrt(s_qq/MC2)*exp(-ybar);
    E_r        = sqrt((1.-x1)*(1.-x2))*MC;
    E0_r       = (sqr(E_r)+sqr(masses[0])-sqr(masses[3]))/(2.*E_r);
    E3_r       = (sqr(E_r)-sqr(masses[0])+sqr(masses[3]))/(2.*E_r);

  } while (E0_r<masses[0] || E3_r<masses[3]);  

  double mt    = sqrt(mt2);
  double pt    = sqrt(mt2-sqr(hadpars.GetConstituents()->Mass(flav)));
  double E1    = mt*cosh(ybar+ystar);
  double pz1   = mt*sinh(ybar+ystar);
  double E2    = mt*cosh(ybar-ystar);
  double pz2   = mt*sinh(ybar-ystar);
  double cosp  = cos(2.*M_PI*ran.Get()), sinp = sqrt(1.-cosp*cosp);
  
  cout<<METHOD<<" : "<<(*cluster)<<endl
      <<"   flav = "<<flav<<" mmax = "<<mmax2
      <<" ("<<MC2<<" "<<masses[0]<<" "<<masses[3]<<"), "
      <<" m_fraction = "<<m_fraction<<" --> "<<endl<<"     "
      <<"mt2, pt, sqq = "<<mt2<<", "<<pt<<", "<<s_qq
      <<", y^* = "<<ystar<<", ybar = "<<ybar<<endl;

  momenta[1]   = Vec4D(E1,pt*cosp,pt*sinp,pz1);
  momenta[2]   = Vec4D(E2,-pt*cosp,-pt*sinp,pz2);
  
  double pz0_r = sqrt(sqr(E0_r)-sqr(masses[0]));
  double gamma = (1.-x1+1.-x2)/(sqrt(4.*(1.-x1)*(1.-x2)));
  double beta  = sqrt(1.-1./sqr(gamma));

  if (ybar<0.) beta *= (-1.);

  double E0    =   E0_r*gamma      - gamma*beta*pz0_r;
  double pz0   = - E0_r*beta*gamma + gamma*pz0_r;
  double E3    =   E3_r*gamma      + gamma*beta*pz0_r;
  double pz3   = - E3_r*beta*gamma - gamma*pz0_r;


  momenta[0]   = Vec4D(E0,0.,0.,pz0);
  momenta[3]   = Vec4D(E3,0.,0.,pz3);

  cout<<"New pair : "<<momenta[1]<<" "<<momenta[1].Abs2()<<"/"
      <<momenta[2]<<" "<<momenta[2].Abs2()<<endl
      <<(momenta[1]+momenta[2]).Abs2()<<" "
      <<(4.*mt2*sqr(cosh(ystar)))<<" "<<(x1*x2*MC2)<<endl;
  
  cout<<"Old pair : "<<momenta[0]<<" "<<momenta[0].Abs2()<<"/"
      <<momenta[3]<<" "<<momenta[3].Abs2()<<endl
      <<(momenta[0]+momenta[3]).Abs2()<<" "<<((1.-x1)*(1.-x2)*MC2)
      <<" from x_{1,2} = "<<x1<<", "<<x2<<", MC2 = "<<MC2<<endl;
   
//    cout<<" for old pair : E_{0,3} = "<<E0_r<<", "<<E3_r<<" and p = "<<pz0_r
//        <<", boost "<<gamma<<","<<beta<<endl;


  cout<<"check this : "<<cluster->Momentum()<<" vs. "
      <<(momenta[0]+momenta[1]+momenta[2]+momenta[3])<<endl
      <<"              axis : "
      <<cluster->Momentum(1)<<"("<<cluster->Momentum(1).Abs2()<<")+"
      <<cluster->Momentum(2)<<"("<<cluster->Momentum(2).Abs2()<<")"<<endl
      <<"                vs.: "
      <<momenta[0]<<"("<<momenta[0].Abs2()<<")+"
      <<momenta[3]<<"("<<momenta[3].Abs2()<<")"<<endl;

  if (m_ana) {
    m_histograms[string("YBar_by_YBarMax")]->Insert(ybar/ybmax);
    m_histograms[string("YBar")]->Insert(ybar);
    m_histograms[string("MT2")]->Insert(mt2);
    m_histograms[string("SQQ")]->Insert(s_qq);
    if (flav==Flavour(kf::u) || flav== Flavour(kf::d)) {
      m_histograms[string("Flavour")]->Insert(0.5);
      m_histograms[string("MT2_{u,d}")]->Insert(mt2);
    }
    else if (flav==Flavour(kf::s)) {
      m_histograms[string("Flavour")]->Insert(1.5);
      m_histograms[string("MT2_s")]->Insert(mt2);
    }
    else if (flav==Flavour(kf::dd_1).Bar() || flav== Flavour(kf::uu_1).Bar() ||
	     flav==Flavour(kf::ud_1).Bar() || flav== Flavour(kf::ud_0).Bar()) {
      m_histograms[string("Flavour")]->Insert(2.5);
      m_histograms[string("MT2_ud")]->Insert(mt2);
    }
    else if (flav==Flavour(kf::sd_1).Bar() || flav== Flavour(kf::su_1).Bar() ||
	     flav==Flavour(kf::sd_0).Bar() || flav== Flavour(kf::su_0).Bar()) {
      m_histograms[string("Flavour")]->Insert(3.5);
      m_histograms[string("MT2_us")]->Insert(mt2);
    }
    else if (flav==Flavour(kf::ss_1).Bar()) {
      m_histograms[string("Flavour")]->Insert(4.5);
      m_histograms[string("MT2_ss")]->Insert(mt2);
    }
  }

  return true;
}

double Schwinger_Uniform::GetYStar(const double mt2,const double mmax2) {
  double ystarmax = 0.5 * log(mmax2/(4.*mt2));
  // allow antiquark forward and backward --> maybe prefer forward only ?
  double ystar, wt, ystarexp = m_ystar_expvalue*ystarmax;
  do {
    ystar = (-1.+2.*ran.Get()) * ystarmax;
    wt    = exp(-sqr((ystar-ystarexp)/m_ystar_sigma));
    //cout<<METHOD<<" "<<ystar<<" from <y> = "<<ystarexp<<" "<<ystarmax<<" "<<wt<<endl;
  } while (wt<ran.Get());

  if (m_ana) {
    m_histograms[string("YStar_by_YStarMax")]->Insert(ystar/ystarmax);
    m_histograms[string("YStar")]->Insert(ystar);
  }
  return ystar;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  QBreakUp Types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool Cluster_Part_QBreakup::BuildKinematics(Cluster * const cluster,Flavour & flav,
					    Vec4D * momenta,double * masses)
{
  double M   = cluster->Mass(0);
  double Q   = SelectQ(M);

  momenta[1] = Q/M*momenta[3];
  momenta[2] = Q/M*momenta[0];
  momenta[0] = (1.-Q/M)*momenta[0];
  momenta[3] = (1.-Q/M)*momenta[3];

  double E1  = sqrt((momenta[0]+momenta[1]).Abs2());
  double E2  = sqrt((momenta[2]+momenta[3]).Abs2());

  flav = Flavour(kf::none);
  bool failure(false);

  do {
    switch (int(m_popper.Pop(flav))) {
    case Return_Value::Success:
    case Return_Value::Warning:
      if (!Veto(cluster,flav)) {
	masses[1] = masses[2] = hadpars.GetConstituents()->Mass(flav);
	if (masses[0]+masses[1]>E1 && masses[2]+masses[3]>E2) {
	  hadpars.AdjustMomenta(4,momenta,masses);
	  return true;
	}
      }
      break;
    case Return_Value::Error: return false;
    }
  } while (!failure);
  return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple Q over M
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Simple_Q_over_M::Simple_Q_over_M() :
  Cluster_Part_QBreakup(0.,int(hadpars.Get(string("FourQ"))))
{ }




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Running Q over M
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Running_Q_over_M::Running_Q_over_M() :
  Cluster_Part_QBreakup(0.,int(hadpars.Get(string("FourQ"))))
{ }

const double Running_Q_over_M::SelectQ(const double M) const
{
  double Q, norm = M*M;
  do { Q = m_Q+ran.Get()*sqrt(m_Q*M); } while (exp(-(Q*Q)/norm)>ran.Get());
  return Q;
}


#include "Cluster_Part.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part() :
  m_4Qmode(0), m_ybar_mode(0), m_ana(true),
  p_popper(hadpars.GetPopper()), 
  m_fraction(sqr(hadpars.Get(string("MassFraction")))),
  m_ystar_expvalue(hadpars.Get(string("<Y>"))),
  m_ystar_sigma(sqr(hadpars.Get(string("Y*_WIDTH"))))
{ 
  if (m_ana) {
    m_histograms[string("YStar_by_YStarMax")] = new Histogram(0,-1.,1.,100);
    m_histograms[string("YStar")]             = new Histogram(0,-5.,5.,100);
    m_histograms[string("YBar_by_YBarMax")]   = new Histogram(0,-1.,1.,100);
    m_histograms[string("YBar")]              = new Histogram(0,-5.,5.,100);
    m_histograms[string("Flavour")]           = new Histogram(0,0.,5.,5);
    m_histograms[string("PT")]                = new Histogram(0,0.,1.5,150);
    m_histograms[string("SQQ")]               = new Histogram(0,0.,5.,100);
  }
}

Cluster_Part::~Cluster_Part()
{
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

bool Cluster_Part::TestDecay(Cluster * const cluster)
{
  //cout<<METHOD<<" before the show "<<endl<<(*cluster)<<endl;
  cluster->BoostInCMSAndRotateOnZ();

  Flavour flav = Flavour(kf::none);
  if (!BuildKinematics(cluster,flav)) return false;

  Cluster * left  = new Cluster(cluster->GetFlav(1),m_momenta[0],flav.Bar(),m_momenta[1]);
  Cluster * right = new Cluster(flav,m_momenta[2],cluster->GetFlav(2),m_momenta[3]);


  if (cluster->GetLeads()==1 || cluster->GetLeads()==3) left->SetLeads(ltp::leadingtrip);
  if (cluster->GetLeads()==2 || cluster->GetLeads()==3) right->SetLeads(ltp::leadinganti);

  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  cluster->RotateAndBoostBack();

  Particle * part;
  part = new Particle(-1,Flavour(kf::cluster),cluster->GetLeft()->Momentum()); 
  part->SetNumber();
  part->SetStatus(part_status::active);
  part->SetInfo('C');
  part->SetFinalMass(cluster->GetLeft()->Mass());
  control::s_AHAparticles++;
  cluster->GetLeft()->SetSelf(part);
  
  part = new Particle(-1,Flavour(kf::cluster),cluster->GetRight()->Momentum()); 
  part->SetNumber();
  part->SetStatus(part_status::active);
  part->SetInfo('C');
  part->SetFinalMass(cluster->GetRight()->Mass());
  control::s_AHAparticles++;
  cluster->GetRight()->SetSelf(part);

  return true;
}


bool Cluster_Part::BuildKinematics(Cluster * const cluster,ATOOLS::Flavour & flav)
{
  double MC    = cluster->Mass(0), MC2 = sqr(MC);
  double m0    = cluster->Mass(1), m02 = sqr(m0);
  double m3    = cluster->Mass(2), m32 = sqr(m3);
  double mmax  = MC - (m0+m3);

  bool diquarks(cluster->GetFlav(1).IsDiQuark()||cluster->GetFlav(2).IsDiQuark());

  flav       = p_popper->SelectFlavour(mmax,diquarks && flav.IsDiQuark());

  if (flav==Flavour(kf::none)) return false;
  double m1    = hadpars.GetConstituents()->Mass(flav), m12 = sqr(m1), m2 = m1, m22 = m12;

  double pmax2  = (sqr(MC2+sqr(m0+m1)-sqr(m2+m3))-4.*sqr(m0+m1)*sqr(m2+m3))/(4.*MC2);
  double ptmax2 = m_fraction*pmax2;
  double pt, ystar, s_qq, ybar, x1, x2, E_r, E0_r(0.), E3_r(0.);
  int counter(0);
  do {
    pt     = p_popper->SelectPT(ptmax2);
    ystar  = GetYStar(sqr(pt)+m12,ptmax2);
    s_qq   = 4.*(sqr(pt)+m12)*sqr(cosh(ystar));
    ybar   = GetYBar(pmax2,s_qq); 
    x1     = sqrt(s_qq/MC2)*exp(ybar);
    x2     = sqrt(s_qq/MC2)*exp(-ybar);
    E_r    = sqrt((1.-x1)*(1.-x2))*MC;
    E0_r   = (sqr(E_r)+m02-m32)/(2.*E_r);
    E3_r   = (sqr(E_r)-m02+m32)/(2.*E_r);

  } while (E0_r<m0 || E3_r<m3);  

  double mt    = sqrt(m12+sqr(pt));
  double E1    = mt*cosh(ybar+ystar);
  double pz1   = mt*sinh(ybar+ystar);
  double E2    = mt*cosh(ybar-ystar);
  double pz2   = mt*sinh(ybar-ystar);
  double cosp  = cos(2.*M_PI*ran.Get()), sinp = sqrt(1.-cosp*cosp);
  
  
  m_momenta[1] = Vec4D(E1,pt*cosp,pt*sinp,pz1);
  m_momenta[2] = Vec4D(E2,-pt*cosp,-pt*sinp,pz2);
  
  double pz0_r = sqrt(sqr(E0_r)-m02);
  double gamma = (1.-x1+1.-x2)/(sqrt(4.*(1.-x1)*(1.-x2)));
  double beta  = sqrt(1.-1./sqr(gamma));

  if (ybar<0.) beta *= (-1.);

  double E0    =   E0_r*gamma      - gamma*beta*pz0_r;
  double pz0   = - E0_r*beta*gamma + gamma*pz0_r;
  double E3    =   E3_r*gamma      + gamma*beta*pz0_r;
  double pz3   = - E3_r*beta*gamma - gamma*pz0_r;


  m_momenta[0] = Vec4D(E0,0.,0.,pz0);
  m_momenta[3] = Vec4D(E3,0.,0.,pz3);


  //   if (cluster->GetFlav(1)==Flavour(kf::c) || cluster->GetFlav(1)==Flavour(kf::b) ||
  //       cluster->GetFlav(2)==Flavour(kf::c).Bar() || cluster->GetFlav(2)==Flavour(kf::b).Bar()) {
  //   cout<<METHOD<<" : "
  //       <<"   flav = "<<flav<<", ptmax = "<<ptmax<<" ("<<MC2<<" "<<m0<<" "<<m3<<"), "
  //       <<" m_fraction = "<<m_fraction<<" --> "<<endl<<"     "
  //       <<"                          pt, sqq = "
  //       <<pt<<", "<<s_qq<<", y^* = "<<ystar<<", ybar = "<<ybar<<endl;
  
  //   cout<<"- New pair              : "
  //       <<m_momenta[1]<<" "<<m_momenta[1].Abs2()<<"/"
  //       <<m_momenta[2]<<" "<<m_momenta[2].Abs2()<<endl
  //       <<"                          "
  //       <<"(p1+p2)^2 = "<<(m_momenta[1]+m_momenta[2]).Abs2()<<" = s_qq = "
  //       <<s_qq<<" vs. "<<(x1*x2*MC2)<<endl;
  
  //   cout<<"- New vecs for old pair : "
  //       <<m_momenta[0]<<" "<<m_momenta[0].Abs2()<<"/"
  //       <<m_momenta[3]<<" "<<m_momenta[3].Abs2()<<endl
  //       <<"                          "
  //       <<(m_momenta[0]+m_momenta[3]).Abs2()<<" "<<((1.-x1)*(1.-x2)*MC2)
  //       <<" from x_{1,2} = "<<x1<<", "<<x2<<", MC2 = "<<MC2<<endl;
  
  //   cout<<" for old pair : E_{0,3} = "<<E0_r<<", "<<E3_r<<" and p = "<<pz0_r
  //       <<", boost "<<gamma<<","<<beta<<endl;
  
  
  //   cout<<"check this : "<<cluster->Momentum()<<" vs. "
  //       <<(m_momenta[0]+m_momenta[1]+m_momenta[2]+m_momenta[3])<<endl
  //       <<"              axis : "
  //       <<cluster->Momentum(1)<<"("<<cluster->Momentum(1).Abs2()<<")+"
  //       <<cluster->Momentum(2)<<"("<<cluster->Momentum(2).Abs2()<<")"<<endl
  //       <<"                vs.: "
  //       <<m_momenta[0]<<"("<<m_momenta[0].Abs2()<<")+"
  //       <<m_momenta[3]<<"("<<m_momenta[3].Abs2()<<")"<<endl;

  //   }

  if (m_ana) {
    m_histograms[string("PT")]->Insert(pt);
    m_histograms[string("SQQ")]->Insert(s_qq);
    if (flav==Flavour(kf::u) || flav== Flavour(kf::d)) {
      m_histograms[string("Flavour")]->Insert(0.5);
    }
    else if (flav==Flavour(kf::s)) {
      m_histograms[string("Flavour")]->Insert(1.5);
    }
    else if (flav==Flavour(kf::dd_1).Bar() || flav== Flavour(kf::uu_1).Bar() ||
	     flav==Flavour(kf::ud_1).Bar() || flav== Flavour(kf::ud_0).Bar()) {
      m_histograms[string("Flavour")]->Insert(2.5);
    }
    else if (flav==Flavour(kf::sd_1).Bar() || flav== Flavour(kf::su_1).Bar() ||
	     flav==Flavour(kf::sd_0).Bar() || flav== Flavour(kf::su_0).Bar()) {
      m_histograms[string("Flavour")]->Insert(3.5);
    }
    else if (flav==Flavour(kf::ss_1).Bar()) {
      m_histograms[string("Flavour")]->Insert(4.5);
    }
  }

  return true;
}

bool Cluster_Part::UpdateDecay(Cluster * const cluster,const int mode)
{
  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(cluster->GetLeft()->GetSelf()->FinalMass()); 
  double m22     = sqr(cluster->GetRight()->GetSelf()->FinalMass());

  cluster->BoostInCMSAndRotateOnZ();
  double pt      = cluster->GetLeft()->Momentum().PPerp();
  double E1      = (M2+m12-m22)/(2.*M);
  double pl1     = sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi  = cos(2.*M_PI*ran.Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

  if (!mode&1) cluster->GetLeft()->RescaleMomentum(p1);
  if (!mode&2) cluster->GetRight()->RescaleMomentum(p2);

  cluster->RotateAndBoostBack();
}


double Cluster_Part::GetYStar(const double mt2,const double mmax2) {
  double ystarmax = 0.5 * log(mmax2/(4.*mt2));
  double ystar, ystarexp = m_ystar_expvalue*ystarmax;
  int maxtrials(100);
  do { 
    ystar = (-1.+2.*ran.Get()) * ystarmax;
    maxtrials--;
  } while (exp(-sqr((ystar-ystarexp)/m_ystar_sigma))<ran.Get() && maxtrials>0);

  if (maxtrials<=0) ystar = (-1.+2.*ran.Get()) * ystarmax;

  if (m_ana) {
    m_histograms[string("YStar_by_YStarMax")]->Insert(ystar/ystarmax);
    m_histograms[string("YStar")]->Insert(ystar);
  }
  return ystar;
}

double Cluster_Part::GetYBar(const double pmax2,const double s_qq) {
  double ybar, ybarmax  = 0.5 * log(pmax2/s_qq);
  double e_ybarmax, cosh_ybarmax;
  int maxtrials(100);
  switch (m_ybar_mode) {
  case 2:
    do { 
      ybar = (-1.+2.*ran.Get());
      maxtrials--;
    } while (exp(-sqr(ybar/ybarmax))<ran.Get() && maxtrials>0);
  case 1:
    e_ybarmax    = exp(ybarmax); 
    cosh_ybarmax = cosh(ybarmax);
    do { 
      ybar = log(1.+ran.Get()*(e_ybarmax-1.));
      if (ran.Get()<.5) ybar *= -1.;
      maxtrials--;
    } while (cosh(ybar)<ran.Get()*cosh_ybarmax && maxtrials>0);
    break;
  case 0:
  default:
    maxtrials = -1;
    break;
  }

  // Uniform distribution
  if (maxtrials<=0) ybar = (-1.+2.*ran.Get())*ybarmax;
  if (m_ana) {
    m_histograms[string("YBar_by_YBarMax")]->Insert(ybar/ybarmax);
    m_histograms[string("YBar")]->Insert(ybar);
  }
  return ybar;
}



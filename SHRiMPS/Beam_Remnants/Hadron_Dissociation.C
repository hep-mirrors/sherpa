#include "SHRiMPS/Beam_Remnants/Hadron_Dissociation.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Hadron_Dissociation::Hadron_Dissociation(Continued_PDF *const pdf) :
  p_pdf(pdf), m_bunch(p_pdf->Bunch()), 
  m_ycut(MBpars("originalY")), //-MBpars("deltaY")),
  m_analyse(true)
{ 
  if (m_analyse) {
    m_histomap[string("KT_remn_orig")] = new Histogram(0,0.0,5.0,50);
    m_histomap[string("KT_remn_resc")] = new Histogram(0,0.0,5.0,50);
    m_histomap[string("X_quark")]      = new Histogram(0,0.0,1.0,1000);
    m_histomap[string("X_gluon")]      = new Histogram(0,0.0,1.0,1000);
    m_histomap[string("X_diquark")]    = new Histogram(0,0.0,1.0,1000);
    m_histomap2D[string("X_quark_2D")] = 
      new Histogram_2D(0,0.,25.,25,0.0,1.0,100);
    m_histomap2D[string("X_gluon_2D")] = 
      new Histogram_2D(0,0.,25.,25,0.0,1.0,100);
    m_histomap2D[string("X_diquark_2D")] = 
      new Histogram_2D(0,0.,25.,25,0.0,1.0,100);
  }
}


Hadron_Dissociation::~Hadron_Dissociation() {
  msg_Info()<<"Initial kt's: "<<m_histomap[string("KT_remn_orig")]->Average()
	    <<" and "<<m_histomap[string("KT_remn_resc")]->Average()
	    <<" after rescaling.\n";
  if (!m_histomap.empty()) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator 
	   hit=m_histomap.begin();hit!=m_histomap.end();hit++) {
      histo = hit->second;
      name  = string("Ladder_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
  if (!m_histomap2D.empty()) {
    Histogram_2D * histo;
    string name;
    for (map<string,Histogram_2D *>::iterator 
	   hit=m_histomap2D.begin();hit!=m_histomap2D.end();hit++) {
      histo = hit->second;
      name  = string("Ladder_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}

void Hadron_Dissociation::FixFlavourConstituents() {
  double random(ran->Get());
  if (m_bunch==Flavour(kf_p_plus)) {
    if (random<1./3.) {
      m_quark   = Flavour(kf_d);
      m_remnant = Flavour(kf_uu_1);
    }     
    else if (random<1./2.) {
      m_quark   = Flavour(kf_u);
      m_remnant = Flavour(kf_ud_1);
    }
    else {
      m_quark   = Flavour(kf_u);
      m_remnant = Flavour(kf_ud_0);
    }
  }
  else if (m_bunch==Flavour(kf_p_plus).Bar()) {
    if (random<1./3.) {
      m_quark   = Flavour(kf_d).Bar();
      m_remnant = Flavour(kf_uu_1).Bar();
    }     
    else if (random<1./2.) {
      m_quark   = Flavour(kf_u).Bar();
      m_remnant = Flavour(kf_ud_1).Bar();
    }
    else {
      m_quark   = Flavour(kf_u).Bar();
      m_remnant = Flavour(kf_ud_0).Bar();
    }
  }
  else {
    msg_Error()<<"Error in "<<METHOD<<"(bunch = "<<m_bunch<<"):\n"
	       <<"   No parton dissociation found.  Will exit.\n";
    exit(1);
  }
}

void Hadron_Dissociation::FillParticleList(const int & N) {
  FixFlavourConstituents();
  msg_Tracking()<<METHOD<<"(N="<<N<<"): "<<m_quark<<" & "<<m_remnant<<".\n";

  m_particles.clear();
  Particle * particle;
  Flavour    defgluon(Flavour(kf_gluon));
  Vec4D      defmom(0.,0.,0.,0.);

  particle = new Particle(0,m_quark,defmom,'B');
  particle->SetNumber(0);
  m_particles.push_back(particle);
  for (int i=0;i<N-1;i++) {
    particle = new Particle(0,defgluon,defmom,'B');
    particle->SetNumber(0);
    m_particles.push_back(particle);
  }
  random_shuffle(m_particles.begin(),m_particles.end(),*ran);

  particle = new Particle(0,m_remnant,defmom,'B');
  particle->SetNumber(0);
  m_particles.push_back(particle);
}

bool Hadron_Dissociation::
DefineDissociation(const int & Nladders,const double B, 
		   const double & xcut,const double & eta,
		   Form_Factor * ff)
{
  //msg_Out()<<METHOD<<"("<<Nladders<<", xcut="<<xcut<<", eta="<<eta<<"):\n";
  m_elastic = false;
  FillParticleList(Nladders+1);

  double  xmin(p_pdf->XMin()*double(Nladders+2)),xave(1./double(Nladders+2));
  double  startweight(pow(1.25,Nladders)*pow(xave,-(2+Nladders)*eta));
  
  if (xmin<xcut) {
    int     trials(0);
    double  weight,wt,x,xsum,maxwt(0.);
    Flavour flav;
    while (trials++<pow(10.,ATOOLS::Min(4,Nladders+3))) {
      m_xs.clear();
      weight = startweight;
      xsum   = 0.;
      int idiquark;
      for (int i=0;i<Nladders+2;i++) {
	xsum += x = xcut + (1.-xcut)*ran->Get();
	m_xs.push_back(x);
      }
      for (int i=0;i<Nladders+2;i++) {
	x    = m_xs[i] /= xsum;
	flav = m_particles[i]->Flav();
	if (flav.IsDiQuark()) {
	  if (x/2<p_pdf->XMin()) { weight = 0.; break; }
          weight *= wt = pow(1.+xcut-x,-Nladders);//exp((x-1.)/Nladders); 
	}
        else if (flav.IsQuark()) {
	  if (x<p_pdf->XMin()) { weight = 0.; break; }
	  p_pdf->Calculate(x,0.);
	  weight *= wt = p_pdf->XPDF(flav)/p_pdf->XPDFMax(flav)/x;
	}
	// the extra x in xpdf compensates for the incoming flux.
	if (i!=Nladders+1 && !IsZero(eta)) weight *= pow(x,eta);
      }
      if (weight>maxwt) maxwt=weight;
      if (weight>1.) {
	msg_Tracking()<<"\n";
	msg_Tracking()<<"   Potential Error in "<<METHOD
		      <<"(Nladders="<<Nladders<<"): "
		      <<"weight = "<<weight<<">1 "
		      <<"(start = "<<startweight<<").\n";
      }
      if (weight>ran->Get()) {
	for (size_t i=0;i<m_particles.size();i++) {
	  m_particles[i]->SetMomentum(m_xs[i] * m_bunchmom);
	  if (m_analyse) {
	    if (m_particles[i]->Flav().IsGluon()){
	      m_histomap["X_gluon"]->Insert(m_xs[i]);
	      m_histomap2D["X_gluon_2D"]->Insert(Nladders,m_xs[i]);
	    }
	    else if (m_particles[i]->Flav().IsQuark()){
	      m_histomap["X_quark"]->Insert(m_xs[i]);
	      m_histomap2D["X_quark_2D"]->Insert(Nladders,m_xs[i]);
	    }
	    else if (m_particles[i]->Flav().IsDiQuark()){
	      m_histomap["X_diquark"]->Insert(m_xs[i]);
	      m_histomap2D["X_diquark_2D"]->Insert(Nladders,m_xs[i]);
	    }
	  }
	}
	DefineTransverseMomenta(ff);
	return true;
      }
    }
    msg_Tracking()<<"\n";
    msg_Error()<<METHOD<<": After "<<trials<<" trials no dissociation for "
    	       <<Nladders<<" ladders, maxwt = "<<maxwt<<".\n";
  }
  return false;
}

void Hadron_Dissociation::Reshuffle(const size_t & N) {
  Particle * help(m_particles[N]);
  double helpx(m_xs[N]);
  Vec4D helpvec(m_qtvecs[N]);
  size_t swap(N);
  if (N==0) swap++; else swap--;
  m_particles[N]    = m_particles[swap];
  m_particles[swap] = help;
  m_xs[N]           = m_xs[swap];
  m_xs[swap]        = helpx;
  m_qtvecs[N]       = m_qtvecs[swap];
  m_qtvecs[swap]    = helpvec;
}

void Hadron_Dissociation::DefineTransverseMomenta(Form_Factor * ff) {
  double QT2min(0.),QT2cut(25.),E,QT2max,KL(0.),qt2,qt,phi;
  Vec4D sumvec(0.,0.,0.,0.);
  vector<double> qts, qls;

  m_qtvecs.clear();
  for (size_t i=0;i<m_xs.size();i++) {
    E      = m_xs[i]*m_bunchmom[0];
    QT2max = Min(QT2cut,E/cosh(m_ycut));
    do { qt2 = ff->SelectQT2(QT2max,QT2min);
    } while (qt2>QT2max || qt2<QT2min);
    qt = sqrt(qt2);
    qts.push_back(qt);
    qls.push_back(sqrt(E*E-qt2));
    if (m_analyse) {
      m_histomap["KT_remn_orig"]->Insert(qt);
    }
    phi     = ran->Get()*2.*M_PI;
    m_qtvecs.push_back(qt*Vec4D(0.,cos(phi),sin(phi),0.));
    sumvec += m_qtvecs.back(); 
    KL     += qls.back();
  }
  for (size_t i=0;i<m_xs.size();i++) {
    m_qtvecs[i] -= qls[i]/KL*sumvec;
    if (m_analyse) {
      qt = m_qtvecs[i].PPerp();
      m_histomap["KT_remn_resc"]->Insert(qt);
    }
  }
}

void Hadron_Dissociation::FillBeamBlob() {
  p_blob->SetType(btp::Beam);
  p_blob->SetTypeSpec("Shrimps");
  p_blob->SetStatus(blob_status::inactive);

  if (!m_elastic) {
    for (size_t i=0;i<m_particles.size();i++) {
      p_blob->AddToOutParticles(m_particles[i]);
    }
  }
  //msg_Out()<<METHOD<<":\n"<<(*p_blob)<<"\n";
}

void Hadron_Dissociation::AddParticlesToBlob(ATOOLS::Blob * blob,int beam) {
  for (size_t i=0;i<m_particles.size();i++) {
    m_particles[i]->SetBeam(beam);
    blob->AddToInParticles(m_particles[i]);
  }
}

void Hadron_Dissociation::PrintParticles() const {
  msg_Out()<<METHOD<<"("<<m_particles.size()<<" particles for "
	   <<m_bunch<<"):\n";
  for (size_t i=0;i<m_particles.size();i++) 
    msg_Out()<<(*m_particles[i])<<"\n";
}

bool Hadron_Dissociation::
MustReplaceColour(const size_t & pos,const size_t & c1,const size_t & c2) {
  for (size_t i=0;i<m_particles.size();i++) {
    if (m_particles[i]->GetFlow(pos)==c1) {
      m_particles[i]->SetFlow(pos,c2);
      return true;
    }
  }
  return false;
}

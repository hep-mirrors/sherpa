#include "BEAM/Spectra/EPA_Spectra_Plotter.H"
#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/EPA_FF.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Special_Functions.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

EPA_Spectra_Plotter::~EPA_Spectra_Plotter() {
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();hit++) {
    histo = hit->second;
    name  = m_dir+string("/")+hit->first+string(".dat");
    msg_Out()<<"Writing "<<name<<"\n";
    histo->Output(name);
    delete histo;
  }
  m_histograms.clear();
  if (p_xaxis) delete p_xaxis;
  if (p_baxis) delete p_baxis;
}

void EPA_Spectra_Plotter::operator()(const size_t disc) {
  switch (disc) {
  case 99:
    TestWoodSaxon();
    TestGauss();
    TestDipole();
    TestPoint();
    break;
  case  5: TestWoodSaxon(); break;
  case  4: TestGauss();     break;
  case  3: TestDipole();    break;
  case  2: TestPoint();     break;
  case  1: TestBesselFunctions(); break;
  default: break;
  }
}

void EPA_Spectra_Plotter::IonFFs() {
  // Still have to fill this ....
}

void EPA_Spectra_Plotter::Init() {
  m_nxbins  = 10000;
  m_nbbins  = 100000;
  if (m_beam.IsIon()) {
    m_xmin  = 1./double(m_nxbins);
    m_xmax  = 1.+m_xmin;
    m_minx  = 1.e-5/m_beam.GetAtomicNumber();
    m_maxx  = 1./m_beam.GetAtomicNumber();
    m_xvals = { 0.001, 0.01, 0.1 };
  }
  else {
    m_xmin  = 1./double(m_nxbins);
    m_xmax  = 1.+m_xmin;
    m_minx  = 1.e-5;
    m_maxx  = 1.;
    m_xvals = { 0.001, 0.01, 0.1 };
  }
  m_xsize   = (m_xmax-m_xmin)/double(m_nxbins);
  m_R       = m_beam.Radius();
  m_bmin    = 0.;
  m_bmax    = 10000.*m_R+m_bmin;
  m_bsize   = (m_bmax-m_bmin)/double(m_nbbins);
  m_nxaxis  = 50;
  m_nbaxis  = 100;
  m_minb    = 1.e-3*m_R*rpa->hBar_c();
  m_maxb    = 1.e6*m_R*rpa->hBar_c();
  p_xaxis   = new axis(m_nxaxis,m_minx,m_maxx,axis_mode::log);
  p_baxis   = new axis(m_nbaxis,m_minb,m_maxb,axis_mode::log);
}

void EPA_Spectra_Plotter::TestPoint() {
  m_beam = Flavour(kf_p_plus);
  Init();
  FillAnalytic(EPA_ff_type::point);
  FillNumerical(EPA_ff_type::point);
  FillNumerical(EPA_ff_type::point,16.);
}

void EPA_Spectra_Plotter::TestDipole() {
  m_beam = Flavour(kf_p_plus);
  Init();
  FillAnalytic(EPA_ff_type::dipole);
  FillNumerical(EPA_ff_type::dipole);
}

void EPA_Spectra_Plotter::TestGauss() {
  m_beam = Flavour(kf_p_plus);
  Init();
  FillAnalytic(EPA_ff_type::Gauss);
  FillNumerical(EPA_ff_type::Gauss);
  m_beam = Flavour(kf_lead208);
  Init();
  FillAnalytic(EPA_ff_type::Gauss);
  FillNumerical(EPA_ff_type::Gauss);
  m_beam = Flavour(kf_calcium40);
  Init();
  FillAnalytic(EPA_ff_type::Gauss);
  FillNumerical(EPA_ff_type::Gauss);
}

void EPA_Spectra_Plotter::TestWoodSaxon() {
  m_beam = Flavour(kf_lead208);
  Init();
  FillNumerical(EPA_ff_type::WoodSaxon);
  m_beam = Flavour(kf_calcium40);
  Init();
  FillNumerical(EPA_ff_type::WoodSaxon);
}


void EPA_Spectra_Plotter::
FillAnalytic(const enum EPA_ff_type & type,const double & Q2max) {
  EPA_FF_Base * ff;
  string ffname;
  switch (type) {
  case EPA_ff_type::WoodSaxon:
    return;
  case EPA_ff_type::dipole:
    ff     = new EPA_Dipole(m_beam, 0);
    ffname = string("dipole");
    break;
  case EPA_ff_type::Gauss:
    ff     = new EPA_Gauss(m_beam, 0);
    ffname = string("gauss");
    break;
  case EPA_ff_type::point:
  default:
    ff     = new EPA_Point(m_beam, 0);
    ff->SetApprox(0);
    ffname = string("point");
    break;
  }
  ff->SetAnalytic(1);
  double pt2max = (Q2max<0. ? sqr(1./m_R) : Q2max);
  ff->SetPT2Max(pt2max);
  string tag           = m_beam.IDName()+string("_")+ffname+string("_analytic");
  if (Q2max>0.) tag   += string("_Q2_")+ToString(pt2max);
  string Nx            = string("Nx_")+tag;
  string Nx_red        = string("Nx_red_")+tag;
  string Ratiox        = string("Ratiox_")+tag;
  m_histograms[Nx]     = new Histogram(0,m_xmin,m_xmax,m_nxbins);
  if (type==EPA_ff_type::point) {
    m_histograms[Nx_red] = new Histogram(0,m_xmin,m_xmax,m_nxbins);
    m_histograms[Ratiox] = new Histogram(0,m_xmin,m_xmax,m_nxbins);
  }
  string Nxb[m_xvals.size()];
  for (size_t i=0;i<m_xvals.size();i++) {
    Nxb[i] = string("Nxb_"+tag+string("_x_")+ToString(m_xvals[i]));
    if (type==EPA_ff_type::point) {
      m_histograms[Nxb[i]] = new Histogram(0,m_bmin,m_bmax,m_nbbins);
    }
  }
  msg_Out()<<METHOD<<" fills histograms with analytic results for "<<m_beam<<":\n"
	   <<"   x in ["<<m_xmin<<", "<<m_xmax<<"], b in ["<<m_bmin<<", "<<m_bmax<<"] "
	   <<"from R = "<<m_R<<".\n";
  for (int xi=0;xi<m_nxbins;xi++) {
    double x   = m_xmin + (double(xi)+0.5)*m_xsize;
    double arg = x / (m_beam.IsIon() ? m_beam.GetAtomicNumber() : 1.);
    double n   = ff->N(arg);
    msg_Out()<<"    N("<<arg<<" --> "<<x<<") = "<<n<<"\n";
    m_histograms[Nx]->Insert(x,n);
    if (type==EPA_ff_type::point) {
      double nred = ff->ReducedN(arg), ratio = nred/n;
      m_histograms[Nx_red]->Insert(x,nred);
      m_histograms[Ratiox]->Insert(x,ratio);
    }
  }
  if (type==EPA_ff_type::point) {
    for (int bi=0;bi<m_nbbins;bi++) {
      double b = m_bmin + (double(bi)+0.5)*m_bsize;
      for (size_t j=0;j<m_xvals.size();j++) {
	double x   = m_xvals[j];
	double arg = x / (m_beam.IsIon() ? m_beam.GetAtomicNumber() : 1.);
	m_histograms[Nxb[j]]->Insert(b,ff->N(arg,b));
      }
    }
  }
  delete ff;
}

void EPA_Spectra_Plotter::
FillNumerical(const EPA_ff_type & type,const double & Q2max) {
  size_t maxapp = 0;
  switch (type) {
  case EPA_ff_type::WoodSaxon: maxapp = 1; break;
  case EPA_ff_type::dipole:    maxapp = 2; break;
  case EPA_ff_type::Gauss:     maxapp = 2; break;
  case EPA_ff_type::point:
  default:                     maxapp = 3; break;
  }
  if (maxapp==0) return;
  EPA_FF_Base * ff[maxapp];
  string partname = m_beam.IDName(), ffname, appname[maxapp], tag[maxapp];
  string Nx[maxapp], Nredx[maxapp], Ratiox[maxapp];
  string Nxb[maxapp][m_xvals.size()];
  for (size_t i=0;i<maxapp;i++) {
    switch (type) {
    case EPA_ff_type::WoodSaxon:
      ff[i] = new EPA_WoodSaxon(m_beam, 0);
      ffname = string("woodsaxon");
      break;
    case EPA_ff_type::dipole:
      ff[i] = new EPA_Dipole(m_beam, 0);
      ffname = string("dipole");
      break;
    case EPA_ff_type::Gauss:
      ff[i] = new EPA_Gauss(m_beam, 0);
      ffname = string("gauss");
      break;
    case EPA_ff_type::point:
    default:
      ff[i] = new EPA_Point(m_beam, 0);
      ffname = string("point");
      break;
    }
    msg_Out()<<"=== "<<METHOD<<"(i = "<<i<<", type = "<<int(type)<<")\n";
    double pt2max = (Q2max<0. ? sqr(1./m_R) : Q2max);
    ff[i]->SetPT2Max(pt2max);
    ff[i]->SetApprox(i);
    ff[i]->SetAnalytic(0);
    //ff[i]->FillTables();
    appname[i] = string("approx_"+ToString(i));
    tag[i]     = partname+string("_")+ffname+string("_")+appname[i];
    if (Q2max>0.) tag[i] += string("_Q2_")+ToString(pt2max);
    Nx[i]                   = string("Nx_")+tag[i];
    Nredx[i]                = string("Nx_red_")+tag[i];
    Ratiox[i]               = string("Ratiox_")+tag[i];
    m_histograms[Nx[i]]     = new Histogram(0,m_xmin,m_xmax,m_nxbins);
    m_histograms[Nredx[i]]  = new Histogram(0,m_xmin,m_xmax,m_nxbins);
    m_histograms[Ratiox[i]] = new Histogram(0,m_xmin,m_xmax,m_nxbins);
    for (size_t j=0;j<m_xvals.size();j++) {
      Nxb[i][j] = string("Nxb_"+tag[i]+"_x_"+ToString(m_xvals[j]));
      m_histograms[Nxb[i][j]] = new Histogram(0,m_bmin,m_bmax,m_nbbins);
    }
  }
  msg_Out()<<METHOD<<" fills histograms with numerical results for "<<m_beam<<"\n"
	   <<"   Output files are named with "<<tag[0]<<"\n";
  for (int bi=0;bi<m_nbbins;bi++) {
    double b = m_bmin + (double(bi)+0.5)*m_bsize;
    for (size_t i=0;i<maxapp;i++) {
      for (size_t j=0;j<m_xvals.size();j++) {
	double x   = m_xvals[j];
	double arg = x / (m_beam.IsIon() ? m_beam.GetAtomicNumber() : 1.);
	m_histograms[Nxb[i][j]]->Insert(b,ff[i]->N(arg,b));
      }
    }
  }
  for (int xi=0;xi<m_nxbins;xi++) {
    double x   = m_xmin + (double(xi)+0.5)*m_xsize;
    double arg = x / (m_beam.IsIon() ? m_beam.GetAtomicNumber() : 1.);
    for (size_t i=0;i<maxapp;i++) {
      double n = ff[i]->N(arg), nred = ff[i]->ReducedN(arg);
      m_histograms[Nx[i]]->Insert(x,n);
      m_histograms[Nredx[i]]->Insert(x,nred);
      m_histograms[Ratiox[i]]->Insert(x,nred/n);
    }
  }
  for (size_t i=0;i<maxapp;i++) delete ff[i];
}

void EPA_Spectra_Plotter::TestBesselFunctions() {
  size_t nbins = 10000;
  double xmin = 0., xmax = 100.+xmin;
  double bsize = (xmax-xmin)/double(nbins);
  msg_Out()<<METHOD<<"(x in "<<"["<<xmin<<", "<<xmax<<"] in "<<nbins<<" bins).\n";
  m_histograms["BesselK0"] = new Histogram(0,xmin,xmax,nbins);
  m_histograms["BesselK1"] = new Histogram(0,xmin,xmax,nbins);
  for (int i=0;i<nbins;i++) {
    double x  = xmin+(double(i)+0.5)*bsize;
    double K0 = SF.Kn(0,x);
    double K1 = SF.Kn(1,x);
    msg_Out()<<x<<" : K0 = "<<K0<<", K1 = "<<K1<<"\n";
    m_histograms["BesselK0"]->Insert(x,K0);
    m_histograms["BesselK1"]->Insert(x,K1);
  }
}

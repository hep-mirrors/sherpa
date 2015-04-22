#include "AHADIC++/Tools/Splitter_Base.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

std::ostream & AHADIC::
operator<<(std::ostream & s, const PoppedPair & pop) {
  s<<" *** ("<<pop.m_flav<<", m = "<<sqrt(pop.m_mpop2)<<") "
   <<"[x = "<<pop.m_x<<", y = "<<pop.m_y<<", z= "<<pop.m_z<<"] "
   <<"[kt = "<<sqrt(pop.m_kt2)<<", mqq = "<<sqrt(pop.m_sqq)<<"]\n"
   <<"     "<<pop.m_outmom[0]<<" ("<<sqrt(Max(0.,pop.m_outmom[0].Abs2()))<<") "
   <<pop.m_outmom[1]<<" ("<<sqrt(Max(0.,pop.m_outmom[1].Abs2()))<<")\n";
  return s;
}

Splitter_Base::Splitter_Base() :
  p_as((Strong_Coupling*)s_model->GetScalarFunction("strong_cpl")),
  m_pt2max(sqr(hadpars->Get(std::string("ptmax")))),
  m_pt02(hadpars->Get(std::string("pt02"))),
  p_split(NULL), p_spect(NULL), m_sumx(0.), m_sumy(0.), m_ana(false)
{
  Init();
  if (m_ana) InitAnalysis();
}

Splitter_Base::~Splitter_Base() {
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) 
    delete fdit->second;
  m_options.clear();

  if (m_ana) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = m_anapath+string("_")+hit->first+std::string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}


void Splitter_Base::Init() {
  double norm(0.);
  Constituents * constituents(hadpars->GetConstituents());
  m_mmin = constituents->MinMass();
  m_mmax = constituents->MaxMass();
  m_mmin2 = ATOOLS::sqr(m_mmin);
  m_mmax2 = ATOOLS::sqr(m_mmax);
  for (FlavCCMap_Iterator fdit=constituents->CCMap.begin();
       fdit!=constituents->CCMap.end();fdit++) {
    if (constituents->TotWeight(fdit->first)>norm)
      norm = constituents->TotWeight(fdit->first);
  }
  DecaySpecs * decspec;
  for (FlavCCMap_Iterator fdit=constituents->CCMap.begin();
       fdit!=constituents->CCMap.end();fdit++) {
    if (!fdit->first.IsAnti()) {
      decspec = new DecaySpecs;
      decspec->popweight = constituents->TotWeight(fdit->first)/norm;
      decspec->massmin   = constituents->Mass(fdit->first);
      m_options.insert(make_pair(fdit->first,decspec));
    }
  }
  if (m_options.empty()) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   No decay channels found for gluons, will abort the run.\n"
	       <<"   Please contact the Sherpa group for assistance.\n";
    exit(0);
  }
  m_sumwt = CalculateSumWT();
}

void Splitter_Base::Reset() {
  while (!m_popped.empty()) {
    delete m_popped.front(); 
    m_popped.pop_front();
  }
  m_sumx = m_sumy = 0.;
}

void Splitter_Base::DefineTags() {
  Flavour flav(p_spect->m_flav);
  m_anti = !((flav.IsQuark() && !flav.IsAnti()) ||
	     (flav.IsDiQuark() && flav.IsAnti()));
  m_leadspect = (p_spect->m_info=='L' || p_spect->m_info=='B');
  m_leadsplit = (p_split->m_info=='L' || p_split->m_info=='B');
  m_isbeam    = (p_split->m_info=='B' || p_spect->m_info=='B');
}

void Splitter_Base::ConstructTrafos() {
  Vec4D cms(p_split->m_mom+p_spect->m_mom);
  m_boost = Poincare(cms);
  m_boost.Boost(p_split->m_mom); 
  m_boost.Boost(p_spect->m_mom);
  m_rotat = Poincare(p_split->m_mom,Vec4D(1.,0.,0.,1.));
  m_rotat.Rotate(p_split->m_mom); 
  m_rotat.Rotate(p_spect->m_mom);
}

void Splitter_Base::UndoTrafos() {
  m_rotat.RotateBack(p_split->m_mom); 
  m_rotat.RotateBack(p_spect->m_mom);
  m_boost.BoostBack(p_split->m_mom); 
  m_boost.BoostBack(p_spect->m_mom);  
}


bool Splitter_Base::ConstructLightC() {
  m_LC.m_smandel = (p_split->m_mom+p_spect->m_mom).Abs2();
  m_LC.m_E       = sqrt(m_LC.m_smandel)/2.;
  m_LC.m_msplit  = p_split->m_flav.HadMass();
  m_LC.m_msplit2 = sqr(m_LC.m_msplit);
  m_LC.m_mspect  = p_spect->m_flav.HadMass();
  m_LC.m_mspect2 = sqr(m_LC.m_mspect);
  if (!AlphaBeta(m_LC.m_smandel,m_LC.m_alpha,m_LC.m_beta)) return false;
  m_LC.m_pA      = m_LC.m_E*Vec4D(1,0.,0.,1);
  m_LC.m_pB      = m_LC.m_E*Vec4D(1,0.,0.,-1);
  Vec4D splitmom((1.-m_LC.m_alpha)*m_LC.m_pA+m_LC.m_beta*m_LC.m_pB);
  Vec4D spectmom(m_LC.m_alpha*m_LC.m_pA+(1.-m_LC.m_beta)*m_LC.m_pB);
  return true;
}

bool Splitter_Base::AlphaBeta(const double & s,double & alpha,double & beta) {
  alpha = beta = 0.;
  if (m_LC.m_msplit>1.e-6 && m_LC.m_mspect>1.e-6) {
    if (sqr(s-m_LC.m_mspect2-m_LC.m_msplit2)<4.*m_LC.m_mspect2*m_LC.m_msplit2) 
      return false;
    double lambda(Lambda(s,m_LC.m_mspect2,m_LC.m_msplit2));
    alpha = (s+m_LC.m_mspect2-m_LC.m_msplit2)/(2.*s)-lambda;
    beta  = m_LC.m_msplit2/(s*(1.-alpha));
  }
  else if (m_LC.m_msplit>1.e-6) {
    alpha = 0.;
    beta  = m_LC.m_msplit2/s;
  }
  else if (m_LC.m_mspect>1.e-6) {
    alpha = m_LC.m_mspect2/s;
    beta  = 0.;
  }
  return (alpha<=1. && beta<=1. && alpha>=-1.e-12 && beta>=-1.e-12);
}

double Splitter_Base::SelectY(const double & ymin,const double & ymax,
			      const double & eta,const double & offset) {
  double y, wt(1.), etap(1.-eta);
  bool logdist(dabs(etap)<=1.e-3);
  double ylow(ymin+offset), yup(ymax+offset);
  do {
    y  = logdist?
      ylow * pow(yup/ylow,ran->Get()):
      pow(pow(ylow,etap)+ran->Get()*(pow(yup,etap)-pow(ylow,etap)),1./etap);
  } while (wt<ran->Get()); // was pow(1.-y,2) before
  return y-offset;
}

double Splitter_Base::SelectZ(const double & delta,const bool & lead) {
  double zmin(0.5*(1.-sqrt(1.-delta))), zmax(0.5*(1.+sqrt(1.-delta))), z;
  do {
    z = zmin+ran->Get()*sqrt(1.-delta);
  } while (1.-2.*z*(1.-z) < ran->Get()); 
  // no structure: flat
  // splitting function: 1.-2.*z*(1.-z)
  // anti splitting[expo] pow(4*z*(1-z),expo)
  return z;
}

bool Splitter_Base::SelectFlavour(const double & sqq,const bool & vetodi) {
  Flavour flav(kf_none);
  double mmax(sqrt(sqq/4.)), m2, sumwt;
  long int calls(0);
  while (calls<100) {
    sumwt = CalculateSumWT(mmax,vetodi);
    if (sumwt>1.e-12) {
      flav = Flavour(kf_none);
      calls++;
      sumwt *= ran->Get();
      for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
	if (fdit->second->popweight>0. && fdit->second->massmin<mmax &&
	    !(vetodi && fdit->first.IsDiQuark())) 
	  sumwt -= fdit->second->popweight;
	if (sumwt<0.) {
	  flav  = fdit->first;
	  m2    = sqr(fdit->second->massmin);
	  break;
	}
      }
      if (PoppedMassPossible(m2)) break;
      mmax = sqrt(m2);
    }
    else {
      for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
	if (vetodi && fdit->first.IsDiQuark()) continue;
      }    
      calls = 100;
    }
  }
  if (calls>=100) {
    //msg_Out()<<"      --> "<<METHOD<<" failed, too many calls.\n"
    //	     <<"      --> sumwt = "<<sumwt<<" for minimal mass = "
    //	     <<sqrt(sqq/4.)<<".\n";
    return false;
  }
  m_popped.back()->m_flav  = flav.IsDiQuark()?flav.Bar():flav; 
  m_popped.back()->m_mpop2 = 
    sqr(hadpars->GetConstituents()->Mass(m_popped.back()->m_flav)); 
  return true;
}

double Splitter_Base::CalculateSumWT(const double & mmax,const bool & vetodi) 
{
  double sumwt(0.), wt;
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
    if (vetodi && fdit->first.IsDiQuark()) continue;
    if (fdit->second->popweight>0. && fdit->second->massmin<0.9999999*mmax) {
      wt = fdit->second->popweight; 
      sumwt += wt;
    }   
  } 
  return sumwt;
}

void Splitter_Base::InitAnalysis() {
  m_histograms[string("y")]   = new Histogram(0,0.,1.,100);
  m_histograms[string("x")]   = new Histogram(0,0.,1.,100);
  m_histograms[string("z")]   = new Histogram(0,0.,1.,100);
  m_histograms[string("kt")]  = new Histogram(0,0.,10.,100);
  m_histograms[string("mqq")] = new Histogram(0,0.,10.,100);
  m_pop = m_popu = m_popd = m_pops = m_popud0 = m_popsu0 = m_popsd0 =
    m_popuu1 = m_popud1 = m_popdd1 = m_popsu1 = m_popsd1 = m_popss1 = 0;
}

void Splitter_Base::Analysis() {
  m_histograms[string("y")]->Insert(m_popped.back()->m_y);
  m_histograms[string("x")]->Insert(m_popped.back()->m_x);
  m_histograms[string("z")]->Insert(m_popped.back()->m_z);
  m_histograms[string("kt")]->Insert(sqrt(Max(0.,m_popped.back()->m_kt2)));
  m_histograms[string("mqq")]->Insert(sqrt(Max(0.,m_popped.back()->m_sqq)));
  m_pop++;
  switch (m_popped.back()->m_flav.Kfcode()) {
  case kf_u:    m_popu++;   break;
  case kf_d:    m_popd++;   break;
  case kf_s:    m_pops++;   break;
  case kf_ud_0: m_popud0++; break;
  case kf_su_0: m_popsu0++; break;
  case kf_sd_0: m_popsd0++; break;
  case kf_uu_1: m_popuu1++; break;
  case kf_ud_1: m_popud1++; break;
  case kf_dd_1: m_popdd1++; break;
  case kf_su_1: m_popsu1++; break;
  case kf_sd_1: m_popsd1++; break;
  case kf_ss_1: m_popss1++; break;
  default: break;
  }
}


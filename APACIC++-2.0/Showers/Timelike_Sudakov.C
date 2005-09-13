#include "Timelike_Sudakov.H"
#include "Run_Parameter.H"
#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Timelike_Kinematics.H"
#include "Sudakov_Tools.H"
#include "Knot.H"
#include "MathTools.H"

#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;
using namespace std;

//-----------------------------------------------------------------------
//-------------------- Constructors -------------------------------------
//----------------------------------------------------------------------- 

Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics * kin,MODEL::Model_Base * model,
				   const double pt2min) :
  p_tools(new Sudakov_Tools(model)), p_kin(kin), m_pt2min(pt2min), m_pt2max(sqr(rpa.gen.Ecms())),
  m_t0(4.*m_pt2min)
{ }

Timelike_Sudakov::~Timelike_Sudakov()
{
  if (p_tools) delete p_tools;
}

void Timelike_Sudakov::Init()
{
  // for static couplings, set first argument to 0.
  p_tools->CalculateMaxCouplings(m_cpl_scheme,m_pt2min,m_pt2max);

  // --- initialise QCD (& QED) splitting functions ---
  for (int i=1;i<17;++i) {
    if (i==7) i=11;
    Flavour fl = Flavour(kf::code(i));
    if (fl.IsOn()) {
      if (fl.Strong()) {
	Add(new q_qg(fl,p_tools));
	Add(new q_qg(fl.Bar(),p_tools));
	if (fl.PSMass()<100.) Add(new g_qq(fl,p_tools));
      }
      if (!(fl.Charge()==0) && (m_direct_photons)) {
	Add(new f_fp(fl,p_tools));
	Add(new f_fp(fl.Bar(),p_tools));
	if (fl.PSMass()<100.) Add(new p_ff(fl,p_tools));
      }
    }
  }
  Add(new g_gg(p_tools));
  
  PrintStat();
}

//-----------------------------------------------------------------------
//-------------------- This Sudakov is a splitting_Group, hence ... -----
//----------------------------------------------------------------------- 

void Timelike_Sudakov::Add(Splitting_Function * spl) 
{
  for (SplFunIter iter(m_group);iter();++iter) {
    if (iter()->GetFlA()==spl->GetFlA()) {
      iter()->Add(spl);
      return ;
    }
  }
  m_group.Append(new Splitting_Group(spl));
  p_selected=spl;
}

double Timelike_Sudakov::CrudeInt(double zmin, double zmax) 
{
  SplFunIter iter(m_group);
  for (;iter();++iter)
    if (iter()->GetFlA()==m_inflav) { p_selected=iter(); break;}
  if (!iter()) return m_lastint = -1.;
  return m_lastint=p_selected->CrudeInt(zmin,zmax);
}        

void Timelike_Sudakov::SelectOne() { p_selected->SelectOne(); }


//-----------------------------------------------------------------------
//-------------------- Dicing the next branch ---------------------------
//----------------------------------------------------------------------- 

bool Timelike_Sudakov::Dice(Knot * mother, Knot * granny) 
{
  //PRINT_INFO("Mother = \n"<<(*mother)<<"\n");
  m_inflav = mother->part->Flav(); 
  m_t      = mother->t;                  
  m_E2     = mother->E2;

  if ((m_t-m_t0)<rpa.gen.Accu()) return 0;

  double z0;

  while (m_t>m_t0) {
    z0  = Max(0.5*(1.-sqrt(1.-m_t0/m_t)),1.e-6);
    CrudeInt(z0,1.-z0);
    //std::cout<<"TL::Dice : Check this out "<<m_t<<"("<<m_t0<<") -> "<<z0<<" -> "<<m_lastint<<std::endl;
    if (m_lastint<0.) return 0;

    if (m_mass_scheme >= 2) m_t -= mother->tout;
    ProduceT(m_t,mother->tout);
    if (m_t<m_t0) return 0;
    if (m_mass_scheme >= 2) m_t += mother->tout;
    // determine estimate for energy
    if (granny) m_E2 = 0.25*sqr(m_t+granny->E2)/granny->E2;
    if (m_t<m_E2) {
      SelectOne();
      m_z   = GetZ();
      //std::cout<<"   TimelikeSudakov : "<<m_t<<" "<<m_z<<" "
      //	       <<m_inflav<<" -> "<<GetFlB()<<" "<<GetFlC()<<std::endl;
      if (!Veto(mother,m_t,m_E2)) {
	if (m_azimuthal_correlation) m_phi = p_selected->GetPhi(m_z);
	                        else m_phi = UniformPhi();
	mother->z      = m_z;
	mother->t      = m_t;
	mother->phi    = m_phi;
	//std::cout<<"Accept emission {"<<m_t<<", "<<m_z<<"} code = "<<m_last_veto<<" for "
	//	 <<mother->kn_no<<"; "<<m_inflav<<std::endl;
	return 1;
      }
      //else std::cout<<"Veto emission {"<<m_t<<", "<<m_z<<"} code = "<<m_last_veto<<std::endl;
    }
  }
  return 0; 
}

//-----------------------------------------------------------------------
//-------------------- Methods for dicing -------------------------------
//----------------------------------------------------------------------- 

void Timelike_Sudakov::ProduceT(double ta, double m2) 
{
  m_t *= exp( 2.*M_PI*log(ran.Get())/m_lastint );
}


bool Timelike_Sudakov::Veto(Knot * mo, double t, double E2) 
{
  // branch okay and mother timelike
  m_last_veto=0;
  // kincheck for first timelike from spacelike
  if (mo->tmax!=0 && m_t>mo->tmax)    return true;
  if (E2<t)                           return true;

  //  sum m1 + m2 < sqrt(ta)
  m_tb  = sqr(GetFlB().PSMass());
  m_tc  = sqr(GetFlC().PSMass());
  m_last_veto=1;
  if (t<sqr(m_tb+m_tc))               return true;
  double z(p_kin->CalcZShift(m_z,t,m_tb,m_tc));
  if (z<0)                            return true;

  m_pt2 = p_kin->CalcKt2(z,m_E2,m_t,m_tb,m_tc);
  if (m_pt2<m_pt2min)                 return true;

  // timelike daughters
  m_last_veto=2;
  double wb(z*z*E2), wc((1.-z)*(1.-z)*E2);
  if ((m_tb>wb) || (m_tc>wc))         return true;

  // z-range and splitting function
  m_last_veto=3;
  if (MassVeto(t,E2))                 return true;

  // 2. alphaS
    m_last_veto=4;
  if (CplVeto())                      return true;

  // 3. angular ordering
  m_last_veto=5;
  if (OrderingVeto(mo,t,E2,z))        return true; 

  // 5. ME
  m_last_veto=7;
  if (MEVeto(mo))                     return true;

  m_last_veto=-1;
  return false;
}

bool Timelike_Sudakov::MassVeto(double t, double E2) 
{
  if (!p_kin->CheckZRange(m_z,E2,t,m_tb,m_tc))                             return true;
  if (GetWeight(m_z,m_pt2,(m_mass_scheme==1||m_mass_scheme==3))<ran.Get()) return true;

  if ((m_width_scheme>0) && (sqr(m_inflav.Width())>0.)) {
    if (m_width_scheme==1 && m_pt2<sqr(m_inflav.Width()))                  return true;
    else if (m_pt2/(m_pt2+sqr(m_inflav.Width()))<ran.Get())                return true;
  }
  return false;
}

bool Timelike_Sudakov::CplVeto() 
{
  switch (m_cpl_scheme) {
  case 0 :
    return false;
  case 2 : 
    return (GetCoupling(0.25*m_t)/GetCoupling()>ran.Get()) ? false : true;   
  default : 
    return (GetCoupling(m_pt2)/GetCoupling()>ran.Get())     ? false : true;   
  }
}

bool Timelike_Sudakov::OrderingVeto(Knot * mo,double t, double E2, double z) 
{
  if (!m_inflav.Strong()) return false;
  switch (m_ordering_scheme) {
  case 0 :                
    return false;
    break;
  case 2 :
    if (m_pt2<mo->maxpt2) return false;
    break;
  case 1 :
  default :
    double th_est(p_kin->CalculateAngle(t,E2,z));
    //std::cout<<"In OrderingVeto for "<<(*mo);
    //if (mo->prev) {
    //  double thprev(p_kin->CalculateAngle(mo->prev->t,mo->prev->E2,mo->prev->z));
    //  std::cout<<"   Prev: "<<thprev<<":"<<(*mo->prev)<<std::endl;
    //}
    if (th_est<mo->thcrit || 3.14<mo->thcrit) {
      mo->costh = th_est;
      return false;
    }
  }
  return true;
}


bool Timelike_Sudakov::MEVeto(Knot * mo) 
{
  if (!m_inflav.Strong()) return false;

  Knot * gr = mo->prev;
  if (gr->t < 0) return false;

  if ((m_MEcorr_scheme == 0) || (!gr))              return false;
  if ((m_MEcorr_scheme == 2) && (gr->prev))         return false;
  if ((m_MEcorr_scheme == 1) && (m_pt2<mo->maxpt2)) return false;
  if (!(m_inflav.IsQuark()))                        return false;

  // determine which is the current twig of the tree:
  Knot * twig = mo;
  while (gr->prev) {
    twig = gr;
    gr   = gr->prev;
  }

  /*
  bool isleft;
  if (twig==gr->left) isleft=1;
                 else isleft=0;
  */

  /*
     if "first" branch perform ME - Correction for gluon radiation
     (t',z') *      
            / \      x_i = 2 E_i / sqrt(t')
     (t,z) *   \
          / \   \
         1   3   2
  */

  double mass123 = gr->t;
  double mass13  = m_t;
  if (m_ordering_scheme == 0) {
    mass123 /= gr->z*(1.-gr->z);
    mass13  /= m_z*(1.-m_z);
  }
  double x2 = 1. - mass13/mass123;
  double x1 = m_z*(2.- x2);           // quark or antiquark !!! 
  double x3 = 2. - x1 -x2;
  double ds_ps =  (1.-x1)/x3 * ( 1. + sqr(x1/(2.-x2)))
                 +(1.-x2)/x3 * ( 1. + sqr(x2/(2.-x1)));
  double ds_me = sqr(x1) + sqr(x2);
  double ratio = ds_me/ds_ps;

  return (ratio<ran.Get()) ? true : false;
}






const Simple_Polarisation_Info & Timelike_Sudakov::GetPolB()
{
  if (m_azimuthal_correlation) m_pol_b = p_selected->GetPolB(m_z,m_phi);
  return m_pol_b;
}

const Simple_Polarisation_Info & Timelike_Sudakov::GetPolC()
{
  if (m_azimuthal_correlation) m_pol_c = p_selected->GetPolC(m_z,m_phi,m_pol_b.Angle());
  return m_pol_c;
}






// ------------------------------------------------------------------------
// --- Ande memorial this part is for consistency checks only -------------
// ------------------------------------------------------------------------

struct Spl_Data {
  int masses;
  double ta;
  std::vector<double> values;
};

void Timelike_Sudakov::CheckSplittings() 
{
  const int dsize=7;
  std::vector<Spl_Data> data(dsize);
  data[0].masses=0;   data[0].ta=0;           // z values
  data[1].masses=0;   data[1].ta=sqr(8.);    // weights
  data[2].masses=1;   data[2].ta=sqr(8.);
  data[3].masses=0;   data[3].ta=sqr(12.);    // weights
  data[4].masses=1;   data[4].ta=sqr(12.);
  data[5].masses=0;   data[5].ta=sqr(50.);    // weights
  data[6].masses=1;   data[6].ta=sqr(50.);
  for (double z=0;z<=1.;z+=10.e-2) {
    data[0].values.push_back(z);
    for (int j=1;j<dsize;++j) {
      data[j].values.push_back(0.);
    }
  }  

  for (SplFunIter iter(m_group);iter();++iter) {
    Splitting_Function * sf = iter();
    sf->PrintStat();
    // calculation
    for (size_t i=1; i<data.size();++i) {
      for (size_t j=0;j<data[0].values.size();++j) {
	double z=data[0].values[j];
	int masses= data[i].masses;
	double ta = data[i].ta;
	m_tb  = sqr(sf->GetFlB().PSMass());
	m_tc  = sqr(sf->GetFlC().PSMass());
	data[i].values[j]=sf->GetWeight(m_z,m_pt2,masses);
      }
    }
    // output 
  }
  THROW(normal_exit,"Finished check.");
}

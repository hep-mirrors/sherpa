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

Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics * kin, double pt2min,
				   MODEL::Model_Base * model, Data_Read * dataread) :
  p_kin(kin), m_pt2min(pt2min) 
{
  /*  angular ordering (1=VO+Coherence, 2=VO)  */ 
  m_ordering_scheme = dataread->GetValue<int>("FS_ORDERING",1);       
  /*  coupling argument (0=fix, 1=pt^2, 2=t/4)  */ 
  m_cpl_scheme      = dataread->GetValue<int>("FS_COUPLINGS",1);
  m_pt_scheme       = dataread->GetValue<int>("FS_PT_DEFINITION",2); 


  /* treatment of massive quarks (0=cuts, 1=a la Catani, 2=define t_eff) */
  m_mass_scheme     = dataread->GetValue<int>("FS_MASS_SCHEME",1);   
  /*  (0=no width supression,  1=cut pt2 always > square of width, 
       2=suppressed by pt2/(Gamma^2+pt2)      */
  m_width_scheme    = dataread->GetValue<int>("FS_WIDTH_SCHEME",0);   
    /*   0=constrained z,  1=unconstrained z)  */
  m_zrange_scheme = dataread->GetValue<int>("FS_Z_RANGES",1);      

  /*  (0=none, 1=hardest so far, 2=first)     */
  m_MEcorr_scheme   = dataread->GetValue<int>("FS_ME_CORRECTIONS",0); 
  /*  (1=approximate angles)  */
  m_angle_scheme    = dataread->GetValue<int>("FS_ANGLE_DEF",1);  
  /*  (0=no photons in shower, 1=photons in shower) */
  m_direct_photons  = dataread->GetValue<int>("FS_PHOTONS",0);        

  m_azimuthal_correlation = dataread->GetValue<int>("FS_ANGLE_CORR",0);

  m_pt2max = sqr(rpa.gen.Ecms());
  p_tools  = new Sudakov_Tools(m_cpl_scheme,model,m_pt2min,m_pt2max);
  m_t0     = 4.*m_pt2min;

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

Timelike_Sudakov::~Timelike_Sudakov()
{
  if (p_tools) delete p_tools;
}


//-----------------------------------------------------------------------
//-------------------- Dicing the next branch ---------------------------
//----------------------------------------------------------------------- 

bool Timelike_Sudakov::Dice(Knot * mother, Knot * granny) 
{
  m_last_veto = 0;
  m_inflav = mother->part->Flav(); 
  m_ta     = mother->t;                  
  m_wa     = mother->E2;

  //  double tend = sqr(sqrt(mother->tout+0.25*m_t0) +sqrt(0.25*m_t0));
  double tend = m_t0;

  if ((m_ta-tend)<rpa.gen.Accu()) return 0;

  while (m_ta>tend) {
    double z0 = 0.5 * (1. - sqrt(1.-m_t0/m_ta));
    //     msg_Debugging()<<" z0(2)="<<z0<<std::endl;
    //     msg_Debugging()<<" z0(1)="<<0.5 * (1. - sqrt(1.-m_t0/m_ta))<<std::endl;
    //    std::cout<<" InitRange: "<<z0<<", "<<1.-z0<<std::endl;
    msg_Debugging()<<" ta="<<m_ta<<std::flush;

    if (z0<rpa.gen.Accu()) {
      msg.Error()<<"In Timelike_Sudakov::Dice : z0 out of bounds : "<<z0<<" !"<<std::endl
		 <<"    Set it on 10^-6 and hope for the best."<<std::endl;
      z0 = 1.e-6;
    }
    CrudeInt(z0,1.-z0);

    if (m_mass_scheme >= 2) m_ta -= mother->tout;
    ProduceT(m_ta,mother->tout);

    double ta_save = m_ta;
    if (m_mass_scheme >= 2) ta_save = m_ta + mother->tout;

    int test_mass=0;
    if ((m_mass_scheme>=2) && (m_ta<m_t0))  test_mass=1;
    if ((m_mass_scheme<2)  && (m_ta<tend))  test_mass=1;
    if (test_mass)  return 0;

    // determine estimate for energy
    if (granny) m_wa = 0.25*sqr(m_ta+granny->E2)/granny->E2;
    double wa_save  = m_wa; 
    if (granny) wa_save = 0.25*sqr(ta_save+granny->E2)/granny->E2;
    double pt2_save = -1.;
    if (ta_save<wa_save) {
      SelectOne();
      m_z   = GetZ();
      
      m_tb  = sqr(GetFlB().PSMass());
      m_tc  = sqr(GetFlC().PSMass());
      double z = p_kin->CalcZShift(m_z,ta_save,m_tb,m_tc);
      if (granny) wa_save = 0.25*sqr(ta_save+granny->E2)/granny->E2;
      pt2_save = m_z*(1.-m_z)*ta_save;
      mother->pt2lcm = pt2_save;
      //      m_pt2    = 0.25*Min((1.-m_z)/m_z,m_z/(1.-m_z))*m_ta; 
      //      m_pt2    = 0.25*Min((1.-z)/z,z/(1.-z))*m_ta; 
      //      std::cout<<" pt2="<<m_pt2<<" ("<<m_tb<<","<<m_tc<<")"<<std::endl;
      m_pt2 = 0.25*p_kin->CalcKt2(z,wa_save,ta_save,m_tb,m_tc);
      //      std::cout<<" M pt2="<<m_pt2<<std::endl;

      if (m_pt2>m_pt2min) {
	//      if (pt2_save>m_pt2min) {
	if (!Veto(mother,ta_save,wa_save)) {
	  m_ta = ta_save; 
	  if (m_azimuthal_correlation) {
	    //	    std::cout<<"calling GetPhi from ["<<mother->kn_no<<"]"<<std::endl;
	    m_phi = p_selected->GetPhi(m_z);
	  }
	  else {
	    UniformPhi();
	  }

	  msg_Debugging()<<mother->part->Flav()<<" Select: z="<<m_z<<" t="<<m_ta<<" pt2="<<m_pt2<<std::endl;
	  mother->z      = m_z;
	  mother->t      = m_ta;
	  mother->phi    = m_phi;
	  if (m_inflav.IsQuark()) mother->maxpt2 = m_pt2;
	  else                    mother->maxpt2 = m_pt2max;
	  return 1;
	}
	msg_Debugging()<<" last_veto="<<m_last_veto<<std::endl;
      }
    }
    m_ta = ta_save; 
  }
  return 0; 
}

//-----------------------------------------------------------------------
//-------------------- Methods for dicing -------------------------------
//----------------------------------------------------------------------- 

void Timelike_Sudakov::ProduceT(double ta, double m2) 
{
  if (m_lastint<0.) m_ta  = -1.;
  else {
    m_ta *= exp( 2.*M_PI*log(ran.Get()) / m_lastint );
  }
}


bool Timelike_Sudakov::Veto(Knot * mo, double ta, double wa) 
{  
  m_last_veto=0;
  double z = m_z;
  if (m_zrange_scheme==1) z = p_kin->CalcZShift(m_z,ta,m_tb,m_tc);

  /* --- to be removed in release version !!! ---
  double phat   = 0.5*(sqrt(wa)+sqrt(wa-ta));
  double ztilde = (m_z*phat + ( (1.-m_z)*ta+m_tb-m_tc)/(4.*phat))/sqrt(wa);
  double ztilde2 = 1.- ((1.-m_z)*phat + ( m_z*ta-m_tb+m_tc)/(4.*phat))/sqrt(wa);
  std::cout<<" t="<<ta<<",  m_tb="<<m_tb<<", m_tc="<<m_tc<<std::endl; 
  std::cout<<" z="<<m_z<<", "<<z<<", "<<ztilde<<","<<ztilde2<<std::endl;

  double E=phat+ta/(4.*phat);
  double E1=m_z*phat+((1-m_z)*ta+m_tb-m_tc)/(4.*phat);
  double E2=(1.-m_z)*phat+(m_z*ta-m_tb+m_tc)/(4.*phat);
  std::cout<<" E1="<<E1<<"  E2="<<E2<<"  E="<<E1+E2<<" ("<<E<<", "<<sqrt(wa)<<")"<<std::endl;

  double zdel=lambda/(2.*ta);
  double z0  =(ta+m_tb-m_tc)/(2.*ta);
  double zm  = z0-zdel;
  double zp  = z0+zdel;
  std::cout<<" zlimits=("<<zm<<","<<zp<<")"<<std::endl;
  */

  double wb      = z*z*wa;
  double wc      = (1.-z)*(1.-z)*wa;
  // timelike daughters
  if ((m_tb>wb) || (m_tc>wc))    return 1;
  // kincheck for first timelike from spacelike
  if (mo->tmax!=0 && m_ta>mo->tmax)  return 1;

  // timelike
  if (wa < ta) {   // (3)
    m_last_veto=1;
    return 1;
  }
  //  sum m1 + m2 < sqrt(ta)
  if (ta  < m_tb+m_tc+2.*sqrt(m_tb*m_tc)) {
    m_last_veto=2;
    return 1;
  }
  // 1. masses, z-range and splitting function
 if (MassVeto(ta,wa)) {  // (2)
   m_last_veto=3;
   return 1;
  }
  // 2. alphaS
  if (CplVeto()) {
    m_last_veto=4;
    return 1;
  }
  // 3. angular ordering
  if (AngleVeto(mo,ta,wa,z)) {
    m_last_veto=5;
    return 1;
  }
  // 4. ME
  if (MEVeto(mo))  {
    m_last_veto=6;
    return 1;
  }
  // 5. JetVeto
//    if (JetVeto(ta,wa,z)) {
//      m_last_veto=7;
//      return 1;    
//    }
  return 0;
}

bool Timelike_Sudakov::MassVeto(double ta, double wa) 
{
  double x_p_mom = sqrt(1.-ta/wa) ;   
  double mean_z,delta_z;              

  switch (m_zrange_scheme) {
  case 0 :
    //    std::cout<<" MassVeto : "<<x_p_mom<<" "<<ta<<std::endl;
    mean_z  = 0.5 *( 1. + (m_tb-m_tc)/ta); 
    delta_z = 0.5 * x_p_mom * sqrt( sqr(ta-m_tb-m_tc) - 4.*m_tb*m_tc )/ta;
    break;
  case 2 :
    mean_z  = 0.5 * (1. + (m_tb-m_tc)/ta); 
    delta_z = 0.5 * sqrt( sqr(ta-m_tb-m_tc) - 4.*m_tb*m_tc )/ta;    
    break;
  case 1 : 
  default:
    mean_z  = 0.5;
    delta_z = 0.5*x_p_mom;
  }
  if ((m_z<mean_z - delta_z) || (m_z>mean_z + delta_z)) {
    msg_Debugging()<<" zrange:"<<mean_z - delta_z<<" < "<<m_z<<" < "<<mean_z + delta_z<<std::endl;
    return 1;
  }
  double w1 = GetWeight(m_z,m_pt2,(m_mass_scheme==1 || m_mass_scheme==3));
  if (w1<ran.Get()) return 1;

  if ((m_width_scheme > 0) && (sqr(m_inflav.Width()) > 0.)) {
    if (m_width_scheme==1) {
      if (m_pt2<sqr(m_inflav.Width())) return 1;
    }
    else {
      if (m_pt2/(m_pt2+sqr(m_inflav.Width())) < ran.Get()) return 1;
    }
  }
  return 0;
}

bool Timelike_Sudakov::CplVeto() 
{
  switch (m_cpl_scheme) {
  case 0 :
    return 0;
    break;
  case 2 : 
    return (GetCoupling(0.25*m_ta)/GetCoupling() > ran.Get()) ? 0 : 1;   
    break;
  default : 
    return (GetCoupling(m_pt2)/GetCoupling() > ran.Get()) ? 0 : 1;   
    break;
  }
}

bool Timelike_Sudakov::AngleVeto(Knot * mo,double ta, double wa, double z) 
{
  if (!m_inflav.Strong()) return 0;
  switch (m_ordering_scheme) {
  case 0 : return 0;
  default :
    {
      double thest = 0.;
      switch (m_angle_scheme) {
      default:  
	thest  = sqrt( ta/(z*(1.- z)*wa) );
	mo->costh = thest;
      }  
      msg_Debugging()<<"check ("<<mo->kn_no<<") th="<<thest<<",   thcrit="<<mo->thcrit<<std::endl;
      
      //if (thest < mo->thcrit) return 0;
      //if (thest < mo->thcrit || (IsEqual(mo->thcrit,M_PI)) ) return 0;
      if (thest < mo->thcrit || (3.14<mo->thcrit && mo->thcrit<3.15) ) return 0;
      return 1;
    }
  }
}

bool Timelike_Sudakov::MEVeto(Knot * mo) 
{
  if (!m_inflav.Strong()) return 0;

  Knot * gr = mo->prev;
  if (gr->t < 0) return 0;

  if ((m_MEcorr_scheme == 0) || (!gr))              return 0;
  if ((m_MEcorr_scheme == 2) && (gr->prev))         return 0;
  if ((m_MEcorr_scheme == 1) && (m_pt2<mo->maxpt2)) return 0;
  if (!(m_inflav.IsQuark()))                        return 0;

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
  double mass13  = m_ta;
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

  return (ratio<ran.Get()) ? 1 : 0;
}


bool Timelike_Sudakov::JetVeto(double ta, double wa, double z) 
{
  if (p_kin->JetVeto(-ta,wa,z,m_tb,m_tc)) {
    msg_Debugging()<<" c JetVeto in Timelike_Kinematics::JetVeto "<<std::endl;
  }

  bool  flag=p_kin->JetVeto(m_ta,m_wa,m_z,0.,0.);
  if (flag) {
    msg_Debugging()<<" d JetVeto in Timelike_Kinematics::JetVeto "<<std::endl;
  }
  return flag;
}


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
    if (iter()->GetFlA()==m_inflav) {p_selected=iter();break;}
  if (!iter()) return m_lastint = -1.;
  return m_lastint=p_selected->CrudeInt(zmin,zmax);
}        

//! Selects one specific mode for the splitting (to be called after CrudeInt)
void Timelike_Sudakov::SelectOne() { p_selected->SelectOne(); }


// ---------------------------------------------------
// --- this part is for consticy checks only ---------
// ---------------------------------------------------

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
	double pt2=z*(1.-z)*ta;
	m_tb  = sqr(sf->GetFlB().PSMass());
	m_tc  = sqr(sf->GetFlC().PSMass());
	if (m_pt_scheme == 1) pt2 -= (1.-m_z)*m_tb + m_z*m_tc;
	else if (m_pt_scheme == 2)
	  m_pt2 = 0.25*Min((1-m_z)/m_z,m_z/(1-m_z))*m_ta;

        data[i].values[j]=sf->GetWeight(m_z,m_pt2,masses);
      }
    }
    // output 
  }
  exit(0);
}


// void Timelike_Sudakov::SetJetvetoPt2(const double pt2) 
// { 
//   p_kin->SetJetvetoPt2(pt2); 
// }

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

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

Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics * _kin,double _pt2min,
				   MODEL::Model_Base * _model,Data_Read * _dataread) :
  p_kin(_kin), m_pt2min(_pt2min) 
{
  m_ordering_scheme = _dataread->GetValue<int>("FS ORDERING",1);       /*  (1=VO+Coherence, 2=VO)                  */ 
  m_cpl_scheme      = _dataread->GetValue<int>("FS COUPLINGS",1);      /*  (0=fix, 1=pt^2, 2=t/4)                  */ 
  m_pt_scheme       = _dataread->GetValue<int>("FS PT DEFINITION",1);  /*  (0=> pt^2 = z(1-z)t      
									   for VO
			                                                   1=> z(1-z)t - (1-z)*t_0(b) - z*t_0(c)   */
  m_mass_scheme     = _dataread->GetValue<int>("FS MASS SCHEME",1);    /*  (0=cuts, 1=a la Catani, 2=define t_eff) */
  m_width_scheme    = _dataread->GetValue<int>("FS WIDTH SCHEME",0);   /*  (0=no width supression,
	 		                                                   1=cut pt2 always > square of width
	  		                                                   2=suppressed by pt2/(Gamma^2+pt2)      */
  if (m_mass_scheme==0) 
    m_zrange_scheme = _dataread->GetValue<int>("FS Z RANGES",1);       //  (only for mass_scheme = 0: 
  else                                                                 //   0=constrained z,
    m_zrange_scheme = 1;                                               //   1=unconstrained z)
  m_MEcorr_scheme   = _dataread->GetValue<int>("FS ME CORRECTIONS",0); /*  (0=none, 1=hardest so far, 2=first)     */
  m_angle_scheme    = _dataread->GetValue<int>("FS ANGLE DEF",1);      /*  (1=approximate angles)                  */
  m_direct_photons  = _dataread->GetValue<int>("FS PHOTONS",0);        /*  (0=no photons in shower,         
									   1=photons in shower)                    */
  m_pt2max = sqr(rpa.gen.Ecms());         // this is an obvious choice ....
  p_tools  = new Sudakov_Tools(m_cpl_scheme,_model,m_pt2min,m_pt2max);
  m_t0     = 4.*m_pt2min;
  //      -- initialise QCD splitting functions -- 
  //      -- initialise QED splitting functions -- 

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

  double tend = sqr(sqrt(mother->tout+0.25*m_t0) +sqrt(0.25*m_t0));

  if ((m_ta-tend)<rpa.gen.Accu()) return 0;
  double z0; 

  while (m_ta>tend) {
    //    if (last_veto==0 || last_veto==7) {
    if (m_mass_scheme >= 2) m_ta -= mother->tout;
    if (m_pt_scheme == 2 ) { z0 = 1. / ( 1. + m_ta/m_t0) ;                                            }
    else                   { z0 = 0.5 * (1. - sqrt(1.-m_t0/m_ta)) ; /* condition that 4 pt^2 > t0 */  }
    if (z0<rpa.gen.Accu()) {
      msg.Error()<<"In Timelike_Sudakov::Dice : z0 out of bounds : "<<z0<<" !"<<std::endl
		 <<"    Set it on 10^-6 and hope for the best."<<std::endl;
      z0 = 1.e-6;
    }
    CrudeInt(z0,1.-z0);

    ProduceT();

    int test_mass=0;
    if ((m_mass_scheme>=2) && (m_ta+mother->tout<tend)) test_mass=1;
    if ((m_mass_scheme<2)  && (m_ta<tend))              test_mass=1;
    if (test_mass) return 0;

    // determine estimate for energy
    if (granny) m_wa = 0.25*sqr(m_ta+granny->E2)/granny->E2;
    if (m_ta<m_wa) {
      SelectOne();
      m_z   = GetZ();
      m_pt2 = m_z*(1.-m_z)*m_ta;
      m_tb  = sqr(GetFlB().PSMass());
      m_tc  = sqr(GetFlC().PSMass());
      if (m_pt_scheme == 1)      m_pt2 -= (1.-m_z)*m_tb + m_z*m_tc;
      else if (m_pt_scheme == 2) m_pt2 = 0.25*Min((1.-m_z)/m_z,m_z/(1.-m_z))*m_ta;
      if (m_pt2>m_pt2min) {
	if (!Veto(mother)) {
	  if (m_mass_scheme >= 2) m_ta += mother->tout;
	  UniformPhi();
	  mother->z      = m_z;
	  mother->t      = m_ta;
	  mother->phi    = m_phi;
	  if (m_inflav.IsQuark()) mother->maxpt2 = m_pt2;
	  else                    mother->maxpt2 = m_pt2max;
	  return 1;
	}
      }
    }
  }
  return 0; 
}

//-----------------------------------------------------------------------
//-------------------- Methods for dicing -------------------------------
//----------------------------------------------------------------------- 

void Timelike_Sudakov::ProduceT() 
{
  if (m_lastint<0.) m_ta  = -1.;
               else m_ta *= exp( 2.*M_PI*log(ran.Get()) / m_lastint );
}


bool Timelike_Sudakov::Veto(Knot * mo) 
{  
  m_last_veto=0;
  double wb      = m_z*m_z*m_wa;
  double wc      = (1.-m_z)*(1.-m_z)*m_wa;
  // timelike daughters
  if ((m_tb>wb) || (m_tc>wc)) {
    return 1;
  }
  // kincheck for first timelike from spacelike
  if (mo->tmax!=0 && m_ta>mo->tmax) {
    return 1;
  }
  // timelike
  if (m_wa < m_ta) {
    m_last_veto=1;
    return 1;
  }
  //  sum m1 + m2 < sqrt(m_ta)
  if (m_ta  < m_tb+m_tc+2.*sqrt(m_tb*m_tc)) {
    m_last_veto=2;
    return 1;
  }
  // 1. masses, z-range and splitting function
  if (MassVeto()) {
    m_last_veto=3;
    return 1;
  }
  // 2. alphaS
  if (CplVeto()) {
    m_last_veto=4;
    return 1;
  }
  // 3. angular ordering
  if (AngleVeto(mo)) {
    m_last_veto=5;
    return 1;
  }
  // 4. ME
  if (MEVeto(mo))  {
    m_last_veto=6;
    return 1;
  }
  // 5. JetVeto
  if (JetVeto(mo)) {
    m_last_veto=7;
    return 1;    
  }
  return 0;
}

bool Timelike_Sudakov::MassVeto() 
{
  double x_p_mom = sqrt(1.-m_ta/m_wa) ;   
  double mean_z,delta_z;              

  switch (m_zrange_scheme) {
  case 0 :
    mean_z  = 0.5 *( 1. + (m_tb-m_tc)/m_ta); 
    delta_z = 0.5 * x_p_mom * sqrt( sqr(m_ta-m_tb-m_tc) - 4.*m_tb*m_tc )/m_ta;
    break;
  default : 
    mean_z  = 0.5;
    delta_z = 0.5*x_p_mom;
    break;
  }
  if ((m_z<mean_z - delta_z) || (m_z>mean_z + delta_z)) return 1;

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

bool Timelike_Sudakov::AngleVeto(Knot * mo) 
{
  if (!m_inflav.Strong()) return 0;
  switch (m_ordering_scheme) {
  case 0 : return 0;
  default :
    double thest;
    switch (m_angle_scheme) {
    default:  
      thest  = sqrt( m_ta/(m_z*(1.- m_z)*m_wa) );
    }  
    if (thest < mo->thcrit) return 0;
    return 1;
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
  bool isleft;
  if (twig==gr->left) isleft=1;
                 else isleft=0;

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


bool Timelike_Sudakov::JetVeto(Knot * mo) 
{
  return p_kin->JetVeto(m_ta,m_wa,m_z,0.,0.);
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

double Timelike_Sudakov::CrudeInt(double _zmin, double _zmax) 
{
  SplFunIter iter(m_group);
  for (;iter();++iter)
    if (iter()->GetFlA()==m_inflav) {p_selected=iter();break;}
  if (!iter()) return m_lastint = -1.;
  return m_lastint=p_selected->CrudeInt(_zmin,_zmax);
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


void Timelike_Sudakov::SetJetvetoPt2(const double pt2) 
{ 
  p_kin->SetJetvetoPt2(pt2); 
}

#include "Timelike_Sudakov.H"
#include "Run_Parameter.H"
#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Timelike_Kinematics.H"

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

//-----------------------------------------------------------------------
//-------------------- Constructors -------------------------------------
//----------------------------------------------------------------------- 

Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics * _kin):kin(_kin) {
  ordering_scheme = 1; /*  (1=VO+Coherence, 2=VO);                           */ 
  cpl_scheme      = 1; /*  (0=fix, 1=pt^2, 2=t/4)                            */ 
  pt_scheme       = 1; /*  (0=> pt^2 = z(1-z)t      for VO
			    1=> z(1-z)t - (1-z)*t_0(b) - z*t_0(c)            */ 	  
  mass_scheme     = 1; /*  (0=cuts, 1=a la Catani, 2=define t_eff)           */
  width_scheme    = 0; /* (1) (0=no width supression,
			    1=cut such that pt2 always > square of width
	  		    2=suppressed according to pt2/(Gamma^2+pt2)      */			    
  zrange_scheme   = 1; /*  (only for mass_scheme = 0: 0=constrained z,
		                                      1=unconstrained z)     */ 
  MEcorr_scheme   = 0; /* "1/0" (0=none, 1=hardest so far, 2=first)          */
  angle_scheme    = 1; /*  (1=approximate angles)                            */
  direct_photons  = 0; /*  (0=no photons in shower, 1=photons in shower)     */

  //      -- initialise alphaS and the crude estimates --
  pt2min = rpa.pshower.FinalQ02();      // rpa.Q02() gives minimum pt^2 for last branch.
  pt2max = sqr(rpa.gen.Ecms());         // this is an obvious choice ....
  tools  = new Sudakov_Tools(cpl_scheme,pt2min,pt2max);
  t0     = 4.*pt2min;

  cout<<" t0="<<t0<<endl;

  //      -- initialise QCD splitting functions -- 
  //      -- initialise QED splitting functions -- 
  for (int i=1;i<17;++i) {
    if (i==7) i=11;
    Flavour fl = Flavour(kf::code(i));
    if (fl.IsOn()) {
      if (fl.Strong()) {
	Add(new q_qg(fl,tools));
	Add(new q_qg(fl.Bar(),tools));
	if (fl.PSMass()<100.) Add(new g_qq(fl,tools));
      }
      if (!(fl.Charge()==0) && (direct_photons)) {
	Add(new f_fp(fl,tools));
	Add(new f_fp(fl.Bar(),tools));
	if (fl.PSMass()<100.) Add(new p_ff(fl,tools));
      }
    };
  }
  Add(new g_gg());

  //  PrintStat();
  // Check splitting functions
  //CheckSplittings();
};

//-----------------------------------------------------------------------
//-------------------- Dicing the next branch ---------------------------
//----------------------------------------------------------------------- 

bool Timelike_Sudakov::Dice(Knot * mother, Knot * granny) {
  inflav = mother->part->Flav(); 
  ta     = mother->t;                  
  wa     = mother->E2;

// old:
//    double tend = t0;
//    if (tend < mother->tout) tend = mother->tout;
// new
  double tend = sqr(sqrt(mother->tout+0.25*t0) +sqrt(0.25*t0));


  msg.Debugging()<<"Timelike_Sudakov::Dice (t,E2): "<<ta<<" / "<<wa<<" / for ("<<mother->kn_no
		 <<"), "<<inflav<<std::endl;

  if ((ta-tend)<rpa.gen.Accu()) {
    msg.Debugging()<<"Timelike_Sudakov::Dice : mother can't branch (ta<t_end) : "
		   <<ta<<" < "<<tend<<std::endl
		   <<"      mother = "<<inflav<<", mass = "<<inflav.PSMass()<<", "
		   <<"status = "<<mother->stat<<std::endl;
    if (mother->prev) {
      msg.Debugging()<<"      prev = "<<mother->prev->part->Flav()  
		     <<", stat = "<<mother->prev->stat<<", t = "<<mother->prev->t<<std::endl;
    }
    return 0; 
  }

  double z0; 

  while (ta>tend) {
    /*
    if ((inflav.IsQuark()) && (mother->tout>t0)) 
      z0 = 0.5*((1. + mother->tout/ta) -
		(1. - mother->tout/ta)*sqrt(1.-(ta*t0)/sqr(ta-mother->tout)));
    else 
    */
    z0 = 0.5 * (1. - sqrt(1.-t0/ta)) ; // condition that pt^2 > t0 !!!

    /* simple checks
       double pa     = 0;
       if (wa>ta) pa=sqrt(wa-ta);
       double ea     = sqrt(wa);
    cout<<" ta =" <<ta<<"\t wa ="<<wa<<"\t z0="<<z0<<endl;
    if (wa>ta) pa=sqrt(wa-ta);
    double z0_test = 0.5 *( 1. - pa/ea);
    cout<<" pa =" <<pa<<"\t ea ="<<ea<<"\t z0="<<z0_test<<endl;
    z0=z0_test;
    */
    
    if (z0<rpa.gen.Accu()) {
      msg.Error()<<"In Timelike_Sudakov::Dice : z0 out of bounds : "<<z0<<" !"<<std::endl;
    }
    CrudeInt(z0,1.-z0);

    if (mass_scheme == 2) ta -= mother->tout;
    ProduceT();
    if (mass_scheme == 2) ta += mother->tout;

    if (ta<tend) {
      msg.Debugging()<<"Timelike_Sudakov::No Branch for ("<<mother->kn_no<<"), "<<inflav
		     <<", set on t="<<mother->tout
		     <<"  fl="<<mother->part->Flav()<<"  mfl="<<mother->part->Flav().PSMass()<<std::endl;
      
      return 0;
    }

    // determine estimate for energy
    if (granny) {
      msg.Debugging()<<"Timelike_Sudakov::Dice() change Energy"<<endl
		     <<" wa = "<<wa;
      wa = 0.25*sqr(ta+granny->E2)/granny->E2;
      msg.Debugging()<<" -> "<<wa<<"  (ta="<<ta<<")"<<endl;
    }

    if (ta<wa) {
      SelectOne();
      z   = GetZ();
      pt2 = z*(1.-z)*ta;
      tb  = sqr(GetFlB().PSMass());
      tc  = sqr(GetFlC().PSMass());
      if (pt_scheme == 1) pt2 -= (1.-z)*tb + z*tc;
      
      if (pt2>pt2min) {
	if (!Veto(mother)) {
	  msg.Debugging()<<"Timelike_Sudakov::Dice Branch with t="<<ta<<", z="<<z<<", ("
			 <<GetFlB()<<","<<GetFlC()<<"), tend="<<tend<<std::endl;      
	  UniformPhi();
	  mother->z      = z;
	  mother->t      = ta;
	  mother->phi    = phi;
	  if (inflav.IsQuark()) mother->maxpt2 = pt2;
	  else mother->maxpt2 = pt2max;
	  return 1;
	}
	else {
	  msg.Tracking()<<"Timelike_Sudakov:  Vetoed Branch with t="<<ta<<", z="<<z<<", ("
			<<GetFlB()<<","<<GetFlC()<<"), tend="<<tend<<std::endl;      
	}
      }
    }
  }
  msg.Debugging()<<"Timelike_Sudakov::Banged out of Dice !"<<std::endl;
  return 0; 
}

//-----------------------------------------------------------------------
//-------------------- Methods for dicing -------------------------------
//----------------------------------------------------------------------- 

void Timelike_Sudakov::ProduceT() {
  if (lastint<0.) ta  = -1.;
             else ta *= exp( 2.*M_PI*log(ran.Get()) / lastint );
}


bool Timelike_Sudakov::Veto(Knot * mo) 
{  
  msg.Debugging()<<std::endl<<"      Enter the vetos with E2, t, z, pt2 = "<<wa<<", "
		 <<ta<<", "<<z<<", "<<z*(1.-z)*ta<<std::endl;
  
  // 0. trivial ranges : timelike, enough energy for daughters and physical opening angle ?

  double wb      = z*z*wa;
  double wc      = (1.-z)*(1.-z)*wa;
  // timelike daughters
  if ((tb>wb) || (tc>wc)) {
    msg.Debugging()<<"      Timelike_Sudakov::Veto : Killed by timelike condition:"
		   <<z<<" "<<wa<<std::endl
		   <<"         d1 "<<tb<<", "<<wb<<" for "<<GetFlB()<<std::endl
		   <<"         d2 "<<tc<<", "<<wc<<" for "<<GetFlC()<<std::endl;
    return 1;
  }

  // timelike
  if (wa < ta) return 1;
  //  sum m1 + m2 < sqrt(ta)
  if (ta  < tb+tc+2.*sqrt(tb*tc)) return 1;

  // 1. masses, z-range and splitting function
  if (MassVeto()) {
    msg.Tracking()<<"MassVeto!"<<endl;
    return 1;
  }
  // 2. alphaS
  if (CplVeto()) {
    msg.Tracking()<<"CplVeto!"<<endl;
    return 1;
  }
  // 3. angular ordering
  if (AngleVeto(mo)) {
    msg.Tracking()<<"AngleVeto!"<<endl;
    return 1;
  }
  // 4. ME
  if (MEVeto(mo))  {
    msg.Tracking()<<"ME!"<<endl;
    return 1;
  }

  // 5. JetVeto *AS*
  if (JetVeto(mo)) {
    msg.Tracking()<<"JetVeto!"<<endl;
    return 1;    
  }
  return 0;
}

bool Timelike_Sudakov::MassVeto() 
{

  // *FK*  double newt0  = sqr( sqrt(tb+0.25*t0) + sqrt(tc+0.25*t0));
  double z0 = 0.5 * (1. - sqrt(1.-t0/ta)) ; // condition that pt^2 > t0 !!!
  double zm = z0;
  double zp = 1. - z0;
  // *FK*  if (z<zm || zp<z) return 1;


  double x_p_mom = sqrt(1.-ta/wa) ;   
  double mean_z,delta_z;              

  switch (zrange_scheme) {
  case 0 :
    mean_z  = 0.5 *( 1. + (tb-tc)/ta); 
    delta_z = 0.5 * x_p_mom * sqrt( sqr(ta-tb-tc) - 4.*tb*tc )/ta;
    break;
  default : 
    mean_z  = 0.5;
    delta_z = 0.5*x_p_mom;
    break;
  }
  if ((z<mean_z - delta_z) || (z>mean_z + delta_z)) {
    //    cout<<" zrange "<<(mean_z - delta_z)<<" < "<<z<<" < "<<(mean_z + delta_z)<<endl;
    return 1;
  }

  double w1 = GetWeight(z,pt2,(mass_scheme==1));
  if (w1<ran.Get()) {
    //    cout<<"weight="<<w1<<endl;
    return 1;
  }

  if ((width_scheme > 0) && (sqr(inflav.Width()) > 0.)) {
    if (width_scheme==1) {
      if (pt2<sqr(inflav.Width())) return 1;
    }
    else {
      if (pt2/(pt2+sqr(inflav.Width())) < ran.Get()) return 1;
    }
  }
  //  msg.Debugging()<<"            MassVeto"<<std::endl;  
  return 0;
}

bool Timelike_Sudakov::CplVeto() {
  msg.Debugging()<<"            In CplVeto. ("<<pt2<<")  "<<std::endl;
  switch (cpl_scheme) {
  case 0 :
    return 0;
    break;
  case 2 : 
    return (GetCoupling(0.25*ta)/GetCoupling() > ran.Get()) ? 0 : 1;   
    break;
  default : 
    return (GetCoupling(pt2)/GetCoupling() > ran.Get()) ? 0 : 1;   
    break;
  }
}

bool Timelike_Sudakov::AngleVeto(Knot * mo) {
  msg.Debugging()<<"            In AngleVeto. ("<<ta<<", "<<z<<")  "<<std::endl;
  if (!inflav.Strong()) return 0;

  switch (ordering_scheme) {
  case 0 : return 0;
  default :
    double thest;
    switch (angle_scheme) {
    default:  
      thest  = sqrt( ta/(z*(1.- z)*wa) );
    }  
    if (thest < mo->thcrit) return 0;
    return 1;
  }
}

bool Timelike_Sudakov::MEVeto(Knot * mo) {
  msg.Debugging()<<"            In MEVeto. ("<<ta<<", "<<z<<")  "<<std::endl;
  if (!inflav.Strong()) return 0;

  Knot * gr = mo->prev;
  if (gr->t < 0) return 0;

  if ((MEcorr_scheme == 0) || (!gr))            return 0;
  if ((MEcorr_scheme == 2) && (gr->prev))       return 0;
  if ((MEcorr_scheme == 1) && (pt2<mo->maxpt2)) return 0;

  // Flavours: ME correction only for q->qg and q->qgamma (has to be done)
  if (!(inflav.IsQuark()))                      return 0;

  // determine which is the current twig of the tree:
  Knot * twig = mo;
  while (gr->prev) {
    twig = gr;
    gr   = gr->prev;
  }
  bool isleft;
  if (twig==gr->left) isleft=1;
                 else isleft=0;

  // if "first" branch perform ME - Correction for gluon radiation
  // (t',z') *      
  //        / \      x_i = 2 E_i / sqrt(t')
  // (t,z) *   \
  //      / \   \
  //     1   3   2

  double mass123 = gr->t;
  double mass13  = ta;
  if (ordering_scheme == 0) {
    mass123 /= gr->z*(1.-gr->z);
    mass13  /= z*(1.-z);
  }
  double x2 = 1. - mass13/mass123;
  double x1 = z*(2.- x2);           // quark or antiquark !!! 
  double x3 = 2. - x1 -x2;

  // without global factor: alpha_s/(2 Pi) * C_F * 1 / ((1-x1) (1-x2));
  double ds_ps =  (1.-x1)/x3 * ( 1. + sqr(x1/(2.-x2)))
                 +(1.-x2)/x3 * ( 1. + sqr(x2/(2.-x1)));
  double ds_me = sqr(x1) + sqr(x2);
  double ratio = ds_me/ds_ps;

  return (ratio<ran.Get()) ? 1 : 0;
}


bool Timelike_Sudakov::JetVeto(Knot * mo) {
  return kin->JetVeto(ta,wa,z,0.,0.);
}


struct Spl_Data {
  int masses;
  double ta;
  std::vector<double> values;
};

void Timelike_Sudakov::CheckSplittings() {
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

  for (SplFunIter iter(group);iter();++iter) {
    Splitting_Function * sf = iter();
    sf->PrintStat();
    // calculation
    for (int i=1; i<data.size();++i) {
      for (int j=0;j<data[0].values.size();++j) {
	double z=data[0].values[j];
	int masses= data[i].masses;
	double ta = data[i].ta;
	double pt2=z*(1.-z)*ta;
	tb  = sqr(sf->GetFlB().PSMass());
	tc  = sqr(sf->GetFlC().PSMass());
	//	if (pt_scheme == 1) pt2 -= (1.-z)*tb + z*tc;
        data[i].values[j]=sf->GetWeight(z,pt2,masses);
      }
    }
    // output
    for (int j=0;j<data[0].values.size();++j) {
      for (int i=0; i<data.size();++i) {
	cout<<data[i].values[j]<<" \t";
      }
      cout<<endl;
    }
 
  };
  exit(0);
}




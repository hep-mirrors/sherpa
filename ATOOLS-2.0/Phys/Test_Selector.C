#include "Test_Selector.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Run_Parameter.H"

using namespace ATOOLS;

Test_Selector::Test_Selector(int _nin,int _nout, Flavour * _fl) { 
  m_nin  = _nin;
  m_nout = _nout;
  m_n    = m_nin+m_nout;

  m_fl = new Flavour[m_nin+m_nout];
  for (int i=0; i<m_nin+m_nout; ++i) {
    m_fl[i]=_fl[i];
  }

  // set cuts (calculate from parameters)

  //Lucifer cuts
  /*
  max_cos_l_l=max_cos_l_beam=0.99619469809;  // i.e. 5 Grad
  max_cos_l_p=max_cos_l_g=0.99619469809; 
  max_cos_p_beam=max_cos_g_beam=max_cos_pg_l=0.99619469809;

  min_m2_q_q = 100.;   // i.e. (10 GeV)^2
  min_E_l=min_E_q=10.; // i.e. 10 GeV
  min_E_p=10.;          // i.e.  1 GeV
  min_E_g=10.;         // i.e. 10 GeV
  */

  //RacoonWW cuts
  
  max_cos_l_l=0.99619469809;     // i.e. 5 Grad
  max_cos_l_beam=0.98480775301;  // i.e. 10 Grad
  max_cos_pg_l=0.99619469809;    // i.e. 5 Grad
  max_cos_p_beam=0.99984769516;  // i.e. 1 Grad 
  max_cos_g_beam=0.99984769516;  // i.e. 1 Grad 

  min_m2_q_q = 25.;   // i.e. (5 GeV)^2
  min_E_l=1.;         // i.e. 1 GeV
  min_E_q=3.;         // i.e.  3 GeV
  min_E_p=0.1;        // i.e. 0.1 GeV
  min_E_g=10.;        // i.e. 10 GeV
  
  //Jegerlehner (hep-ph/0109290)
  /*  
  max_cos_l_l= 1.;    // i.e. 0 Grad
  max_cos_l_beam=0.985;          // i.e. 10 Grad
  max_cos_pg_l=0.99619469809;    // i.e. 5 Grad
  max_cos_p_beam=0.985;          // i.e. 10 Grad 
  max_cos_g_beam= 1.;  // i.e. 1 Grad 

  min_m2_q_q = 100.;   // i.e. (10 GeV)^2
  //min_m2_l_l = 0.;   // i.e. (10 GeV)^2
  min_E_l=5.;          // i.e. 5 GeV
  min_E_q=0.;          // i.e.  0 GeV
  min_E_p=1.;           // i.e. 1 GeV
  min_E_g=0.;         // i.e. 10 GeV
  */
  sel_logs.push_back(new Selector_Log("min_E_g")); // 0
  sel_logs.push_back(new Selector_Log("min_E_p")); // 1
  sel_logs.push_back(new Selector_Log("min_E_l")); // 2
  sel_logs.push_back(new Selector_Log("min_E_q")); // 3
  sel_logs.push_back(new Selector_Log("max_cos_l_beam")); // 4
  sel_logs.push_back(new Selector_Log("max_cos_p_beam")); // 5
  sel_logs.push_back(new Selector_Log("max_cos_g_beam")); // 6
  sel_logs.push_back(new Selector_Log("min_m2_q_q")); // 7
  sel_logs.push_back(new Selector_Log("max_cos_l_l")); // 8
  sel_logs.push_back(new Selector_Log("max_cos_pg_l")); // 9
}

bool Test_Selector::Trigger(const Vec4D* vecs) 
{
  cur_vecs=vecs;
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    // minimal Energy of gluon and photon
    if ((m_fl[i].IsGluon())&&(vecs[i][0]<min_E_g)) 
      return sel_logs[0]->Hit();
    if ((m_fl[i].Kfcode()==kf::photon)&&(vecs[i][0]<min_E_p)) 
      return sel_logs[1]->Hit();


    // minimal Energy of charged lepton and quark
    if ((m_fl[i].IntCharge()!=0)&&(m_fl[i].IsLepton())
	&&(vecs[i][0]<min_E_l)) return sel_logs[2]->Hit();
    if ((m_fl[i].IsQuark())&&(vecs[i][0]<min_E_q)) 
      return sel_logs[3]->Hit();
    
    // minimal angle between lepton/quark/photon/gluon <-> beam
    for (int k=0; k<m_nin ; ++k) {
      if ((m_fl[i].IntCharge()!=0)&&(m_fl[i].IsFermion())&&(!m_fl[i].IsQuark())
	  &&(cos_ij(i,k)>max_cos_l_beam)) return sel_logs[4]->Hit();   
      if ((m_fl[i].Kfcode()==kf::photon)&&(cos_ij(i,k)>max_cos_p_beam)) {
	return sel_logs[5]->Hit();   
      }
      if ((m_fl[i].IsGluon())&&(cos_ij(i,k)>max_cos_g_beam)) 
	return sel_logs[6]->Hit();   
    }

    for (int j=i+1; j<m_nin+m_nout; ++j) {
  
      if ((m_fl[j].IsQuark())&&(m_fl[i].IsQuark())) {
	// minmal invariant mass squared between quarks
	if (m2_ij(i,j)<min_m2_q_q) return sel_logs[7]->Hit();   
      } else if (m_fl[i].IsFermion()  && m_fl[j].IsFermion() &&
		 m_fl[i].IntCharge()!=0 && m_fl[j].IntCharge()!=0) {
	// minimal angle between lepton/quark <-> lepton
	if (cos_ij(i,j)>max_cos_l_l) return sel_logs[8]->Hit(); 
      } else {
	// photon/gluon <-> quark/lepton
	if ((m_fl[i].IntCharge()!=0)&&(m_fl[i].IsFermion())&&
	    ((m_fl[j].Kfcode()==kf::photon)||(m_fl[j].IsGluon()))&&
	      (cos_ij(i,j)>max_cos_pg_l)) return sel_logs[9]->Hit();  
	if (((m_fl[i].Kfcode()==kf::photon)||(m_fl[i].IsGluon()))&&
	    (m_fl[j].IntCharge()!=0)&&(m_fl[j].IsFermion())&&
	    (cos_ij(i,j)>max_cos_pg_l)) return sel_logs[9]->Hit();  
      }
    } 
  }
  return 1;
}


void Test_Selector::BuildCuts(Cut_Data* cuts) 
{
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    // minimal Energy of gluon and photon
    if (m_fl[i].IsGluon())            cuts->energymin[i] = min_E_g; 
    if (m_fl[i].Kfcode()==kf::photon) cuts->energymin[i] = min_E_p;

    // minimal Energy of charged lepton and quark
    if ((m_fl[i].IntCharge()!=0)&&(m_fl[i].IsLepton())) cuts->energymin[i] = min_E_l;

    if (m_fl[i].IsQuark()) cuts->energymin[i] = min_E_q; 
    
    // minimal angle between lepton/quark/photon/gluon <-> beam
    for (int k=0; k<m_nin ; ++k) {
      if ((m_fl[i].IntCharge()!=0)&&(m_fl[i].IsFermion())&&(!m_fl[i].IsQuark())) 
	cuts->cosmax[i][k] = cuts->cosmax[k][i] = max_cos_l_beam;

      if (m_fl[i].Kfcode()==kf::photon) 
	cuts->cosmax[i][k] = cuts->cosmax[k][i] = max_cos_p_beam;

      if (m_fl[i].Kfcode()==kf::gluon) 
	cuts->cosmax[i][k] = cuts->cosmax[k][i] = max_cos_g_beam;
    }
    
    for (int j=i+1; j<m_nin+m_nout; ++j) {
  
      if ((m_fl[j].IsQuark())&&(m_fl[i].IsQuark())) {
	// minimal invariant mass squared between quarks
	cuts->scut[i][j] = cuts->scut[j][i] = min_m2_q_q; 
      } else if (m_fl[i].IsFermion()  && m_fl[j].IsFermion() &&
		 m_fl[i].IntCharge()!=0 && m_fl[j].IntCharge()!=0) {
	// minimal angle between lepton/quark <-> lepton
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = max_cos_l_l;
      } else {
	// photon/gluon <-> quark/lepton
	if ((m_fl[i].IntCharge()!=0)&&(m_fl[i].IsFermion())&&
	    ((m_fl[j].Kfcode()==kf::photon)||(m_fl[j].IsGluon()))) 
	  cuts->cosmax[i][j] = cuts->cosmax[j][i] = max_cos_pg_l;
      }
    } 
  }


  //minimal energy = mass
  for (int i=0; i<m_nin+m_nout; ++i) cuts->energymin[i] = Max(cuts->energymin[i],rpa.consts.Mass(m_fl[i],sqr(rpa.gen.Ecms())));

  //scut
  for (int i=0; i<m_nin+m_nout; ++i) {
    for (int j=i+1; j<m_nin+m_nout; ++j) {
      double sc = 
	+sqr(rpa.consts.Mass(m_fl[i],sqr(rpa.gen.Ecms())))+sqr(rpa.consts.Mass(m_fl[j],sqr(rpa.gen.Ecms())))
	+2.*cuts->energymin[i]*cuts->energymin[j]
	-2.*sqrt(dabs(sqr(cuts->energymin[i])-sqr(rpa.consts.Mass(m_fl[i],sqr(rpa.gen.Ecms())))))
  	   *sqrt(dabs(sqr(cuts->energymin[j])-sqr(rpa.consts.Mass(m_fl[j],sqr(rpa.gen.Ecms())))))
	*cuts->cosmax[i][j];
      cuts->scut[i][j] = Max(Max(cuts->scut[i][j],sc),1.e-12*sqr(rpa.gen.Ecms()));
      cuts->scut[j][i] = cuts->scut[i][j];
    }
  } 
}











































#include "Test_Selector.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Run_Parameter.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS; // scut init
using namespace AMATOOLS;

Test_Selector::Test_Selector(int _Nin,int _Nout, Flavour * _flavs) : 
    Nin(_Nin), Nout(_Nout) {

  flavs = new Flavour[_Nin+_Nout];
  for (int i=0; i<Nin+Nout; ++i) {
    flavs[i]=_flavs[i];
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

bool Test_Selector::Trigger(const vec4d* vecs) 
{
  cur_vecs=vecs;
  for (int i=Nin; i<Nin+Nout; ++i) {
    // minimal Energy of gluon and photon
    if ((flavs[i].isgluon())&&(vecs[i][0]<min_E_g)) 
      return sel_logs[0]->Hit();
    if ((flavs[i].kfcode()==kf::photon)&&(vecs[i][0]<min_E_p)) 
      return sel_logs[1]->Hit();


    // minimal Energy of charged lepton and quark
    if ((flavs[i].icharge()!=0)&&(flavs[i].islepton())
	&&(vecs[i][0]<min_E_l)) return sel_logs[2]->Hit();
    if ((flavs[i].isquark())&&(vecs[i][0]<min_E_q)) 
      return sel_logs[3]->Hit();
    
    // minimal angle between lepton/quark/photon/gluon <-> beam
    for (int k=0; k<Nin ; ++k) {
      if ((flavs[i].icharge()!=0)&&(flavs[i].isfermion())&&(!flavs[i].isquark())
	  &&(cos_ij(i,k)>max_cos_l_beam)) return sel_logs[4]->Hit();   
      //if ((flavs[i].kfcode()==kf::photon)) cout<<"p_beam : "<<cos_ij(i,k)<<endl;
      if ((flavs[i].kfcode()==kf::photon)&&(cos_ij(i,k)>max_cos_p_beam)) {
	//cout<<"p_beam (kicked) : "<<cos_ij(i,k)<<endl;
	return sel_logs[5]->Hit();   
      }
      if ((flavs[i].isgluon())&&(cos_ij(i,k)>max_cos_g_beam)) 
	return sel_logs[6]->Hit();   
    }

    for (int j=i+1; j<Nin+Nout; ++j) {
  
      if ((flavs[j].isquark())&&(flavs[i].isquark())) {
	// minmal invariant mass squared between quarks
	if (m2_ij(i,j)<min_m2_q_q) return sel_logs[7]->Hit();   
      } else if (flavs[i].isfermion()  && flavs[j].isfermion() &&
		 flavs[i].icharge()!=0 && flavs[j].icharge()!=0) {
	// minimal angle between lepton/quark <-> lepton
	if (cos_ij(i,j)>max_cos_l_l) return sel_logs[8]->Hit(); 
      } else {
	// photon/gluon <-> quark/lepton
	if ((flavs[i].icharge()!=0)&&(flavs[i].isfermion())&&
	    ((flavs[j].kfcode()==kf::photon)||(flavs[j].isgluon()))&&
	      (cos_ij(i,j)>max_cos_pg_l)) return sel_logs[9]->Hit();  
	if (((flavs[i].kfcode()==kf::photon)||(flavs[i].isgluon()))&&
	    (flavs[j].icharge()!=0)&&(flavs[j].isfermion())&&
	    (cos_ij(i,j)>max_cos_pg_l)) return sel_logs[9]->Hit();  
      }
    } 
  }
  //  cout<<" done "<<endl;
  return 1;
}


void Test_Selector::BuildCuts(Cut_Data* cuts) 
{
  for (int i=Nin; i<Nin+Nout; ++i) {
    // minimal Energy of gluon and photon
    if (flavs[i].isgluon())            cuts->energymin[i] = min_E_g; 
    if (flavs[i].kfcode()==kf::photon) cuts->energymin[i] = min_E_p;

    // minimal Energy of charged lepton and quark
    if ((flavs[i].icharge()!=0)&&(flavs[i].islepton())) cuts->energymin[i] = min_E_l;

    if (flavs[i].isquark()) cuts->energymin[i] = min_E_q; 
    
    // minimal angle between lepton/quark/photon/gluon <-> beam
    for (int k=0; k<Nin ; ++k) {
      if ((flavs[i].icharge()!=0)&&(flavs[i].isfermion())&&(!flavs[i].isquark())) 
	cuts->cosmax[i][k] = cuts->cosmax[k][i] = max_cos_l_beam;

      if (flavs[i].kfcode()==kf::photon) 
	cuts->cosmax[i][k] = cuts->cosmax[k][i] = max_cos_p_beam;

      if (flavs[i].kfcode()==kf::gluon) 
	cuts->cosmax[i][k] = cuts->cosmax[k][i] = max_cos_g_beam;
    }
    
    for (int j=i+1; j<Nin+Nout; ++j) {
  
      if ((flavs[j].isquark())&&(flavs[i].isquark())) {
	// minimal invariant mass squared between quarks
	cuts->scut[i][j] = cuts->scut[j][i] = min_m2_q_q; 
      } else if (flavs[i].isfermion()  && flavs[j].isfermion() &&
		 flavs[i].icharge()!=0 && flavs[j].icharge()!=0) {
	// minimal angle between lepton/quark <-> lepton
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = max_cos_l_l;
      } else {
	// photon/gluon <-> quark/lepton
	if ((flavs[i].icharge()!=0)&&(flavs[i].isfermion())&&
	    ((flavs[j].kfcode()==kf::photon)||(flavs[j].isgluon()))) 
	  cuts->cosmax[i][j] = cuts->cosmax[j][i] = max_cos_pg_l;
      }
    } 
  }

  // used

  /*
  cout<<" in void Test_Selector::BuildCuts(Cut_Data* cuts) :"<<endl;

  for (int i=0;i<Nin+Nout; ++i) {
    for (int j=i+1;j<Nin+Nout; ++j) {
      cout<<" cut["<<j<<"]["<<i<<"]="<<cuts->scut[j][i]<<endl;      
    }
  }
  */

  //minimal energy = mass
  for (int i=0; i<Nin+Nout; ++i) cuts->energymin[i] = Max(cuts->energymin[i],rpa.consts.Mass(flavs[i],sqr(rpa.gen.Ecms())));

  //scut
  for (int i=0; i<Nin+Nout; ++i) {
    for (int j=i+1; j<Nin+Nout; ++j) {
      double sc = 
	+sqr(rpa.consts.Mass(flavs[i],sqr(rpa.gen.Ecms())))+sqr(rpa.consts.Mass(flavs[j],sqr(rpa.gen.Ecms())))
	+2.*cuts->energymin[i]*cuts->energymin[j]
	-2.*sqrt(dabs(sqr(cuts->energymin[i])-sqr(rpa.consts.Mass(flavs[i],sqr(rpa.gen.Ecms())))))
  	   *sqrt(dabs(sqr(cuts->energymin[j])-sqr(rpa.consts.Mass(flavs[j],sqr(rpa.gen.Ecms())))))
	*cuts->cosmax[i][j];
      cuts->scut[i][j] = Max(Max(cuts->scut[i][j],sc),1.e-12*sqr(rpa.gen.Ecms()));
      //cuts->scut[i][j] = 1.e-14*sqr(rpa.gen.Ecms());
      cuts->scut[j][i] = cuts->scut[i][j];
      //      cout<<i<<";"<<j<<" : "<<cuts->scut[i][j]<<endl;
    }
  } 

  /*
  cout<<" in void Test_Selector::BuildCuts(Cut_Data* cuts) :"<<endl;

  for (int i=0;i<Nin+Nout; ++i) {
    for (int j=i+1;j<Nin+Nout; ++j) {
      cout<<" cut["<<j<<"]["<<i<<"]="<<cuts->scut[j][i]<<endl;      
    }
  }
  */

}



/* old
int Test_Selector::Trigger(const vec4d* vecs) {
  cur_vecs=vecs;
  for (int i=Nin; i<Nin+Nout; ++i) {
    // minimal Energy of gluon and photon
    if ((flavs[i].isgluon())&&(vecs[i][0]<min_E_g)) return 0;
    if ((flavs[i].kfcode()==kf::photon)&&(vecs[i][0]<min_E_p)) return 0;


    if (((flavs[i].islepton())||
	(flavs[i].isquark())) &&
	(flavs[i].icharge()!=0)) {  // charged lepton/quark
      // minimal Energy of charged lepton and quark
      if ((flavs[i].islepton())&&(vecs[i][0]<min_E_l)) return 0;
      if ((flavs[i].isquark())&&(vecs[i][0]<min_E_q)) return 0;
    
      // minimal angle between lepton/quark <-> beam
      if (cos_ij(i,0)>max_cos_l_beam) return 0;   


      for (int j=i+1; j<Nin+Nout; ++j) {
	if (((flavs[j].isquark())||(flavs[j].islepton())) &&
	    (flavs[j].icharge()!=0)) {	  
	  
	  if ((flavs[j].isquark())&&(flavs[i].isquark())) {
	    // minmal invariant mass squared between quarks
	    if (m2_ij(i,j)<min_m2_q_q) return 0;   
	  } else {
	    // minimal angle between lepton/quark <-> lepton
	    if (cos_ij(i,j)>max_cos_l_l) return 0; 
	  }
	} 
      }
    }
  }
  return 1;
}
*/








































#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "HADRONS++/PS_Library/Two_Body_PSs.H"
#include "HADRONS++/PS_Library/Three_Body_PSs.H"
#include "HADRONS++/PS_Library/Four_Body_PSs.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"

using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


bool HD_Channel_Selector::DecomposeChannel(string name, ChannelInformation & ci)
{
  ci.name = "noname";
  ci.a=0; ci.b=0; ci.c=0; ci.d=0;
  ci.res1 = "no res";
  ci.res2 = "no res";
  ci.res3 = "no res";
  
  Data_Reader reader("_",";","#");
  vector<string> exploded;
  reader.SetString(name);
  reader.VectorFromString(exploded);
  
  if(exploded.size() < 1) return false;
  
  if(exploded[0]=="Isotropic" || exploded[0]=="Iso2") {
    ci.name = exploded[0];
    ci.nRes = 0;
  }
  else if(exploded[0]=="Dalitz" && exploded.size()==3) {
    ci.name = exploded[0];
    ci.res1 = exploded[1];
    int ab = ToType<int>(exploded[2]);
    ci.b=ab%10; ci.a=ab/10;   // int/int !
    ci.nRes = 1;
  }
  else if(exploded[0]=="TwoResonances" && exploded.size()==5) {
    ci.name = exploded[0];
    ci.res1 = exploded[1]; 
    ci.a    = ToType<int>(exploded[2]);
    ci.res2 = exploded[3]; 
    int bc = ToType<int>(exploded[4]);
    ci.c=bc%10; ci.b=bc/10;   // int/int !
    ci.nRes = 2;
  }
  else if(exploded[0]=="IsotropicSpectator" && exploded.size()==2) {
    ci.name = exploded[0];
    ci.a = ToType<int>(exploded[1]); // spectator index
  }
  
  if( ci.name==string("noname") ) return false;
  return true;
}

Single_Channel * HD_Channel_Selector::GetChannel( 
    int nin, 
    int nout, 
    const Flavour * flavs, 
    string name,
    GeneralModel const & md,
    const ATOOLS::Mass_Selector* ms)
{
  if ( nin>1 || nout<1 ) {
    msg_Error()<<METHOD<<": Error: "<<endl
           <<"   No PS for channel ("<<nin<<" -> "<<nout<<" )"<<endl
           <<"   Return nothing and hope for the best."<<endl;
    return NULL;
  }
  ChannelInformation ci;
  if( DecomposeChannel( name, ci ) ) {
    if (ci.name==string("Isotropic")) {
      if ( nout == 2 ) return new Iso2Channel(flavs);
      if ( nout == 1 ) return new Iso1Channel(flavs);
      return new Rambo(1,nout,flavs,ms);
    }
  }
  if (ci.name==string("Iso2") || nout==2 ) return new Iso2Channel(flavs);
  if (nout==3) {
    if (ci.name==string("Dalitz")) {
      kf_code kfres (kf_rho_770_plus);
      if( ci.res1==string("photon") ) kfres = kf_photon;
      if( ci.res1==string("rho(770)+") ) kfres = kf_rho_770_plus;
      if( ci.res1==string("K*(892)+") ) kfres = kf_K_star_892_plus;
      if( ci.res1==string("rho(1700)+") ) kfres = kf_rho_1700_plus;
      if( ci.res1==string("J/psi(1S)") ) kfres = kf_J_psi_1S;
      if( ci.res1==string("psi(2S)") ) kfres = kf_psi_2S;
      if( ci.res1==string("psi(4040)") ) kfres = kf_psi_4040;
      double width =  Flavour(kfres).Width();
      if( ci.res1==string("W") ) {
        kfres = kf_Wplus;
        width = 2.06;
      }
      SimpleResonanceFlavour res(
          Flavour(kfres).IDName(),
          md("Mass_"+Flavour(kfres).IDName(), Flavour(kfres).HadMass() ),
          md("Width_"+Flavour(kfres).IDName(), width ) );
      return new Dalitz(flavs,res,ci.a,ci.b);
    }
  }
  if (nout==4) {
    if( ci.name==string("TwoResonances") ) {
      SimpleResonanceFlavour res_a( 
          ci.res1, 
          md("Mass_"+ci.res1, Flavour(kf_a_1_1260_plus).HadMass()),
          md("Width_"+ci.res1,Flavour(kf_a_1_1260_plus).Width())); 
      string helpname = ci.res2; // vector resonanance as it appears in md
      SimpleResonanceFlavour res_v( 
          ci.res2,
          md("Mass_"+helpname, Flavour(kf_rho_770_plus).HadMass()),
          md("Width_"+helpname,Flavour(kf_rho_770_plus).Width()) ); 
      return new TwoResonances( flavs, res_a, ci.a, res_v, ci.b, ci.c );
    }
    if( ci.name==string("IsotropicSpectator") ) {
      return new IsotropicSpectator( flavs, ci.a, ms );
    }
  }

  msg_Error()<<METHOD<<": Error: "<<endl
    <<"   No channel for ("<<nin<<" -> "<<nout<<") with name "<<name<<endl
         <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}

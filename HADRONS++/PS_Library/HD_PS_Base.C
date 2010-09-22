#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/PS_Library/Two_Body_PSs.H"
#include "HADRONS++/PS_Library/Three_Body_PSs.H"
#include "HADRONS++/PS_Library/Four_Body_PSs.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/Message.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"

using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

#include "ATOOLS/Org/MyStrStream.H"

////////// class HD_Channel_Selector /////////
 
bool HD_Channel_Selector::DecomposeChannel( string name, ChannelInformation & ci )
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
    GeneralModel const & md )
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
      return new Rambo(1,nout,flavs,true);
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
      if( ci.res1==string("W") ) kfres = kf_Wplus;
      SimpleResonanceFlavour res(
          Flavour(kfres).IDName(),
          md("Mass_"+Flavour(kfres).IDName(), Flavour(kfres).HadMass() ),
          md("Width_"+Flavour(kfres).IDName(), Flavour(kfres).Width() ) );
      return new Dalitz(flavs,res,ci.a,ci.b);
    }
  }
  if (nout==4) {
    if( ci.name==string("TwoResonances") ) {
      SimpleResonanceFlavour res_a( 
          ci.res1, 
          md("Mass_"+ci.res1, Flavour(kf_a_1_1260_plus).HadMass()),
          md("Width_"+ci.res1,Flavour(kf_a_1_1260_plus).Width())); 
      string helpname;                      // name of vector resonanance as it appears in md
      helpname = ci.res2;                   // take name unchanged
      SimpleResonanceFlavour res_v( 
          ci.res2,
          md("Mass_"+helpname, Flavour(kf_rho_770_plus).HadMass()),
          md("Width_"+helpname,Flavour(kf_rho_770_plus).Width()) ); 
      return new TwoResonances( flavs, res_a, ci.a, res_v, ci.b, ci.c );
    }
    if( ci.name==string("IsotropicSpectator") ) {
      return new IsotropicSpectator( flavs, ci.a );
    }
  }

  msg_Error()<<METHOD<<": Error: "<<endl
    <<"   No channel for ("<<nin<<" -> "<<nout<<") with name "<<name<<endl
         <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}

////////// class HD_PS_Base /////////

HD_PS_Base::HD_PS_Base( Hadron_Decay_Channel * hdc ) :
  Multi_Channel("hadron decay channel"),
  p_channelselector(new HD_Channel_Selector), p_hdc(hdc),
  m_res(-1.), m_error(1.), m_max(-1.), m_flux(1./(2.*hdc->Flavours()[0].HadMass()))
{
  nin=1;
  nout=hdc->NOut();
}


HD_PS_Base::~HD_PS_Base()
{
  delete p_channelselector; p_channelselector=NULL;
}

bool HD_PS_Base::IsChannel( string name )
{
  GeneralModel ghost_md;
  Single_Channel * sc = p_channelselector->GetChannel( 1, p_hdc->NOut(),
                                                       p_hdc->Flavours(),
                                                       name, ghost_md );
  if (sc==NULL) return 0;
  delete sc;
  return 1;
}

bool HD_PS_Base::AddChannel(string name,double weight, GeneralModel const & md)
{
  Single_Channel * sc = p_channelselector->GetChannel( 1, p_hdc->NOut(),
                                                       p_hdc->Flavours(),
                                                       name, md );
  if (sc==NULL) return 0;
  sc->SetAlpha(weight);
  Add(sc);                                  // add this to channels in Multi_Channel
  return 1;
}

vector<double> HD_PS_Base::CalculateNormalisedWidth() {
  msg_Info()<<METHOD<<" for "<<p_hdc->Name()<<endl;
  Reset();
  long int iter = Number()*5000*int(pow(2.,int(p_hdc->NOut())-2));
  int maxopt    = Number()*int(pow(2.,2*(int(p_hdc->NOut())-2)));

  long int n=0;
  int      opt=0;
  double   value, oldvalue=0., sum=0., sum2=0., result=1., disc;
  bool     simple=false;

  while(opt<maxopt && result>0. && m_error/result>0.002 ) {
    for (n=1;n<iter+1;n++) {
      value = p_hdc->Differential();
      sum  += value;
      sum2 += ATOOLS::sqr(value);
      AddPoint(value);
      if (value>m_max) {
        m_max = value;
      }
      if (value!=0. && value==oldvalue) { simple = true; break; }
      oldvalue = value;
    }
    opt++;
    Optimize(0.01);

    if (simple) break;          // this way error=0
    n      = opt*iter;
    result = sum/n;
    disc   = sqr(sum/n)/((sum2/n - sqr(sum/n))/(n-1));
    if (disc!=0.0) m_error  = result/sqrt(abs(disc));
    msg_Info()<<"     result (w/o flux): "<<result<<" +/- "<<m_error<<" ("
              <<m_error/result*100.<<" %)"<<endl;
  }
  m_res  = m_flux*sum/n;
  m_error *= m_flux;
  disc   = sqr(m_res)/((sum2*sqr(m_flux)/n - sqr(m_res))/(n-1));
  if (disc!=0.0) m_error  = m_res/sqrt(abs(disc));
  if(abs(m_error)/m_res<1e-6) m_error=0.0;
  msg_Info()<<"     result (incl. flux): "<<m_res<<" +/- "<<m_error<<" ("<<m_error/m_res*100.<<" %)"<<endl;
  // note: the m_max is w/o flux factor

  vector<double> results;
  results.push_back(m_res);
  results.push_back(m_error);
  results.push_back(m_max);
  return results;
}


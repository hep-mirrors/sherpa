#include "HD_PS_Base.H"
#include "Hadron_Decay_Channel.H"
#include "Two_Body_PSs.H"
#include "Three_Body_PSs.H"
#include "Four_Body_PSs.H"
#include "Rambo.H"
#include "Data_Reader.H"
#include "Message.H"
#include "ResonanceFlavour.H"

using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

#include "MyStrStream.H"

////////// class HD_Channel_Selector /////////
 
bool HD_Channel_Selector::DecomposeChannel( string name, ChannelInformation & ci )
{
  char s[name.size()];
  strcpy( s, name.c_str() );
  char delim[] = "_";
  char *result (NULL);
  result = strtok( s, delim );
  int i (0);
  ci.name = "noname";
  ci.a=0; ci.b=0; ci.c=0; ci.d=0;
  ci.res1 = "no res";
  ci.res2 = "no res";
  ci.res3 = "no res";
  while( result != NULL ) {
    if( strcmp(result,"Isotropic" )==0 ||
        strcmp(result,"Iso2")==0 ) {
      ci.name = result;
      ci.nRes = 0;
    }
    if( strcmp(result,"Dalitz")==0 ) {
      ci.name = result;
      result = strtok( NULL, delim ); ci.res1 = result;
      result = strtok( NULL, delim ); int ab = atoi( result );
      ci.b=ab%10; ci.a=ab/10;   // int/int !
      ci.nRes = 1;  
    }
    if( strcmp(result,"TwoResonances")==0 ) {
      ci.name = result;
      result = strtok( NULL, delim ); ci.res1 = result; 
      result = strtok( NULL, delim ); ci.a = atoi( result );
      result = strtok( NULL, delim ); ci.res2 = result; 
      result = strtok( NULL, delim ); int bc = atoi( result );
      ci.c=bc%10; ci.b=bc/10;   // int/int !
      ci.nRes = 2;
    }
    result = strtok( NULL, delim );
  }
  if( ci.name==string("noname") ) return 0;
  return 1;
}

Single_Channel * HD_Channel_Selector::GetChannel( 
    int nin, 
    int nout, 
    const Flavour * flavs, 
    string name,
    GeneralModel const & md )
{
  if (flavs[0].Kfcode() == kf::K ) return NULL;
  if ( nin>1 || nout<2 ) {
    msg.Error()<<METHOD<<": Error: "<<endl
           <<"   No PS for channel ("<<nin<<" -> "<<nout<<" )"<<endl
           <<"   Return nothing and hope for the best."<<endl;
    return NULL;
  }
  ChannelInformation ci;
  if( DecomposeChannel( name, ci ) ) {
    if (ci.name==string("Isotropic")) {
      if ( nout == 2 ) return new Iso2Channel(flavs);
      return new Rambo(1,nout,flavs);
    }
  }
  if (ci.name==string("Iso2") || nout==2 ) return new Iso2Channel(flavs);
  if (nout==3) {
    if (ci.name==string("Dalitz")) {
      kf::code kfres (kf::rho_770_plus);
      if( ci.res1==string("photon") ) kfres = kf::photon;
      if( ci.res1==string("rho(770)+") ) kfres = kf::rho_770_plus;
      if( ci.res1==string("K*(892)+") ) kfres = kf::K_star_892_plus;
      if( ci.res1==string("rho(1700)+") ) kfres = kf::rho_1700_plus;
      if( ci.res1==string("W") ) kfres = kf::W;
      SimpleResonanceFlavour res(
          Flavour(kfres).IDName(),
          md("Mass_"+Flavour(kfres).IDName(), Flavour(kfres).PSMass() ),
          md("Width_"+Flavour(kfres).IDName(), Flavour(kfres).Width() ) );
      return new Dalitz(flavs,res,ci.a,ci.b);
    }
  }
  if (nout==4) {
    if( ci.name==string("TwoResonances") ) {
      SimpleResonanceFlavour res_a( 
          ci.res1, 
          md("Mass_"+ci.res1, Flavour(kf::a_1_1260_plus).PSMass()),
          md("Width_"+ci.res1,Flavour(kf::a_1_1260_plus).Width())); 
      string helpname;                      // name of vector resonanance as it appears in md
      helpname = ci.res2;                   // take name unchanged
      SimpleResonanceFlavour res_v( 
          ci.res2,
          md("Mass_"+helpname, Flavour(kf::rho_770_plus).PSMass()),
          md("Width_"+helpname,Flavour(kf::rho_770_plus).Width()) ); 
      return new TwoResonances( flavs, res_a, ci.a, res_v, ci.b, ci.c );
    }
  }

  msg.Error()<<METHOD<<": Error: "<<endl
    <<"   No channel for ("<<nin<<" -> "<<nout<<") with name "<<name<<endl
         <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}

////////// class HD_PS_Base /////////

HD_PS_Base::HD_PS_Base( Hadron_Decay_Channel * hdc ) :
  Multi_Channel("hadron decay channel"), p_hdc(hdc),
  p_channelselector(new HD_Channel_Selector),
  m_res(-1.), m_error(0.), m_max(-1.), m_flux(1./(2.*hdc->Flavours()[0].Mass()))
{
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
  msg.Info()<<"HD_PS_Base::CalculateNormalisedWidth() for "
    <<p_hdc->ChannelName()<<endl;
  Reset();
  long int iter = Number()*5000*int(pow(2.,int(p_hdc->NOut())-2));
  int maxopt    = Number()*int(pow(2.,2*(int(p_hdc->NOut())-2)));

  long int n;
  int      opt=0;
  double   value, oldvalue=0., sum=0., sum2=0., result=-1., disc;
  bool     simple=false;
  bool isotropic_me = false; // fixme
  if( p_hdc->GetCurrents().size()==0 &&
      p_hdc->GetMEs().size()==1 &&
      p_hdc->GetMEs()[0].second->METype() == "Isotropic" ) isotropic_me = true;

  while(opt<maxopt || (result>0. && m_error/result>0.01) ) {
    for (n=1;n<iter+1;n++) {
      value = p_hdc->Differential(NULL,NULL);
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
    if (disc>0) m_error  = result/sqrt(disc);
    msg.Info()<<"     result (w/o flux): "<<result<<" +/- "<<m_error<<" ("<<m_error/result*100.<<" %)"<<endl;
    if (isotropic_me && m_error/result < 0.01) break;
    if(m_error/result < 0.0007) break;
  }
  m_res  = m_flux*sum/n;
  m_error *= m_flux;
  disc   = sqr(m_res)/((sum2*sqr(m_flux)/n - sqr(m_res))/(n-1));
  if (disc>0) m_error  = m_res/sqrt(disc);
  msg.Info()<<"     result (incl. flux): "<<m_res<<" +/- "<<m_error<<" ("<<m_error/m_res*100.<<" %)"<<endl;
  // note: the m_max is w/o flux factor

  vector<double> results;
  results.push_back(m_res);
  results.push_back(m_error);
  results.push_back(m_max);
  return results;
}



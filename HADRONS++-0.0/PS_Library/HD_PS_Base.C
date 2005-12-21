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
    GeneralModel & md )
{
  if (nin>1 || nout<2) {
    msg.Error()<<"Error in HD_Channel_Selector::GetChannel : "<<endl
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
      ResonanceFlavour res;
      if( ci.res1==string("photon") )    res.Set( kf::photon, 0., 0. );
      if( ci.res1==string("rho(770)") )  
        res.Set( 
            kf::rho_770, 
            md("Mass_Rho_770",  Flavour(kf::rho_770_plus).PSMass()), 
            md("Width_Rho_770", Flavour(kf::rho_770_plus).Width()) );
      return new Dalitz(flavs,res,ci.a,ci.b);
    }
  }
  if (nout==4) {
    if( ci.name==string("TwoResonances") ) {
      ResonanceFlavour res_a( 
          ci.res1, 
          md("Mass_"+ci.res1, Flavour(kf::a_1_1260_plus).PSMass()),
          md("Width_"+ci.res1,Flavour(kf::a_1_1260_plus).Width())); 
      string helpname;                      // name of resonanance as it appears in md
      helpname = ci.res2;                   // take name unchanged
      if( (int)helpname[helpname.size()-1] >= 48 &&
          (int)helpname[helpname.size()-1] <= 57 ) {    // if last char is a number
        helpname.insert( helpname.size()-1, "_" );      // insert _ inbetween
      }
      ResonanceFlavour res_v( 
          ci.res2,
          md("Mass_"+helpname, Flavour(kf::rho_770_plus).PSMass()),
          md("Width_"+helpname,Flavour(kf::rho_770_plus).Width()) ); 
      return new TwoResonances( flavs, res_a, ci.a, res_v, ci.b, ci.c );
    }
  }

  msg.Error()<<"Error in HD_Channel_Selector::GetChannel : "<<endl
    <<"   No channel for ("<<nin<<" -> "<<nout<<") with name "<<name<<endl
         <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}

////////// class HD_PS_Base /////////

HD_PS_Base::HD_PS_Base( 
    Hadron_Decay_Channel * hdc, 
    string _path,
    vector<string> & _pst, 
    bool & mustinit, 
    struct GeneralModel &_locmd,
    bool read_dc ) :
  Multi_Channel("hadron decay channel"), p_hdc(hdc),
  p_channelselector(new HD_Channel_Selector), m_file(string("")),
  m_res(-1.), m_error(0.), m_max(-1.), m_flux(1./(2.*hdc->Flavours()[0].Mass())),
  m_read_dcfile( read_dc ),
  m_path(_path), m_foundPS( false )
{
  if (_pst.size()>2) m_file = _pst[2];      // filename of DC file
  mustinit = Construct(_locmd);             // call Construct to do the rest
  delete p_channelselector;
}


HD_PS_Base::~HD_PS_Base() {}

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

void HD_PS_Base::Initialise() 
{
  CalculateNormalisedWidth();
  WriteOut();
}

bool HD_PS_Base::Construct( GeneralModel & _md )
{
  bool mustinit (true);
  if (m_file!=string("") && m_read_dcfile) {        // in case there is a DC file => read it !
    msg_Tracking()<<"HD_PS_Base::Construct(...) : read "<<m_path<<m_file<<endl;
    vector<vector<string> > helpsvv;
    Data_Reader reader = Data_Reader(string("|"),string(";"),string("!"));
    reader.SetAddCommandLine(false);
    reader.AddComment("#");
    reader.AddComment("//");
    reader.SetInputPath(m_path);
    reader.SetInputFile(m_file);
    reader.SetMatrixType(mtc::transposed);
    if(!reader.MatrixFromFile(helpsvv)) {
      msg.Error()<<"ERROR in HD_PS_Base::Construct(...) :\n"
         <<"   Read in failure "<<m_path<<m_file<<", will abort."<<endl;
      abort();
    }
    
    string name;
    double weight;
    bool skipresult (false);
    vector<string> channels_vector;
    vector<double> ch_weights;
    channels_vector.clear();
    ch_weights.clear();
    for (int i=0;i<helpsvv.size();i++) {
      if ( helpsvv[i][0]==string("Channels") ) {
        i++;                                            // next line
        while (helpsvv[i][0]!=string("}")) {
          weight=1.;
          if (helpsvv[i][0]==string("AlwaysIntegrate")) {
            skipresult = atoi( helpsvv[i][2].c_str() );
          }
          else {
            if (helpsvv[i].size()>1) {                  // if factor is given
              weight=atof(helpsvv[i][1].c_str());
            }
            else weight = 1.;
            if (IsChannel(helpsvv[i][0])) {             // if it is a channel that Sherpa can cope with
              channels_vector.push_back( helpsvv[i][0] );       // save it for later
              ch_weights.push_back( weight );
            }
          }
          i++;
        }
        m_foundPS = 1;                              
        if (channels_vector.size()==0) {                // no channel found
          AddChannel( string("Isotropic"), 1., _md );   //   take Rambo
          skipresult = true;                            //   and don't read Result
          m_foundPS = 0;                
        }
      }
      if ( helpsvv[i][0] == string("Resonances") ) {
        i++;
        while (helpsvv[i][0]!=string("}")) {
          if ( helpsvv[i][1] == string("->") ) {
            if ( helpsvv[i][0] == string("vector1") ) {         // particles into which V1 decays
              _md["vector1_i"] = ToType<double>(reader.Interpreter()->Interprete(helpsvv[i][2]));
              _md["vector1_j"] = ToType<double>(reader.Interpreter()->Interprete(helpsvv[i][3]));
            }
            if ( helpsvv[i][0] == string("vector2") ) {         // particles into which V2 decays
              _md["vector2_i"] = ToType<double>(reader.Interpreter()->Interprete(helpsvv[i][2]));
              _md["vector2_j"] = ToType<double>(reader.Interpreter()->Interprete(helpsvv[i][3]));
            }
            i++;
          }
        }
      }
      if ( helpsvv[i][0] == string("Parameters") ) {

        // in DC file: complex values are to be given in "abs" "phase"
        i++;
        while ( helpsvv[i][0] != string("}") ) {
          if ( helpsvv[i][1] == string("=")) {
            if( helpsvv[i].size() == 3 ) {        // <name> = <real value>
              _md[helpsvv[i][0]] = ToType<double> (
                  reader.Interpreter()->Interprete(helpsvv[i][2]) );
            }
            if( helpsvv[i].size() == 4 ) {        // <name> = <complex value>
              _md[helpsvv[i][0]+string("_abs")] = ToType<double> (
                  reader.Interpreter()->Interprete(helpsvv[i][2]) );
              _md[helpsvv[i][0]+string("_phase")] = ToType<double> (
                  reader.Interpreter()->Interprete(helpsvv[i][3]) );
            }
          }
          if ( helpsvv[i][2] == string("=")) {
            if( helpsvv[i].size() == 4 ) {        // <name> <index> = <real value>
              _md[helpsvv[i][0]+string("_")+helpsvv[i][1]] = ToType<double> (
                  reader.Interpreter()->Interprete(helpsvv[i][3]) );
            }
            if( helpsvv[i].size() == 5 ) {        // <name> <index> = <complex value>
              _md[helpsvv[i][0]+string("_")+helpsvv[i][1]+string("_abs")] 
                = ToType<double> ( reader.Interpreter()->Interprete(helpsvv[i][3]) );
              _md[helpsvv[i][0]+string("_")+helpsvv[i][1]+string("_phase")] 
                = ToType<double> ( reader.Interpreter()->Interprete(helpsvv[i][4]) );
            }
          }
          i++;
        }
      }
      if (helpsvv[i][0]==string("Result") && !skipresult ) {
        i++;
        while (helpsvv[i][0]!=string("}")) {
          m_res   = atof(helpsvv[i][0].c_str());
          m_error = atof(helpsvv[i][1].c_str());
          m_max   = atof(helpsvv[i][2].c_str());
          i++;
        }
        mustinit = false;   // return: no need to integrate
      }
    }
    for (int i=0; i<channels_vector.size(); i++) {
      AddChannel( channels_vector[i], ch_weights[i], _md ); 
    }
  }
  else {                    // in case there is no DC file
    AddChannel( string("Isotropic"), 1., _md );
  }
  return mustinit;          // return: it has to be integrated or not
}

bool HD_PS_Base::AddChannel(string name,double weight,GeneralModel & md) 
{
  Single_Channel * sc = p_channelselector->GetChannel( 1, p_hdc->NOut(), 
                                                       p_hdc->Flavours(),
                                                       name, md );
  if (sc==NULL) return 0;
  sc->SetAlpha(weight);
  Add(sc);                                  // add this to channels in Multi_Channel
  return 1;
}

void HD_PS_Base::CalculateNormalisedWidth() {
  msg.Out()<<"HD_PS_Base::CalculateNormalisedWidth() for "
    <<p_hdc->ChannelName()<<endl;
  Reset();
  long int iter = Number()*5000*int(pow(2.,int(p_hdc->NOut())-2));
  int maxopt    = Number()*int(pow(2.,2*(int(p_hdc->NOut())-2)));
  
  long int n;
  int      opt=0;
  double   value, oldvalue=0., sum=0., sum2=0., result=-1., disc;
  bool     maxincrease, simple=false;
  bool     isotropic_me = (p_hdc->GetME()->METype() == "Isotropic" )? 1 : 0;
  while(opt<maxopt || (result>0. && m_error/result>0.01) ) {
    maxincrease = false;
    for (n=1;n<iter+1;n++) {
      value = p_hdc->Differential();
      sum  += value;
      sum2 += ATOOLS::sqr(value);
      AddPoint(value);
      if (value>m_max) { 
        m_max = value; maxincrease = true; 
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
  } 
  m_res  = m_flux*sum/n;
  m_max *= m_flux;
  m_error *= m_flux;
  disc   = sqr(m_res)/((sum2*sqr(m_flux)/n - sqr(m_res))/(n-1));
  if (disc>0) m_error  = m_res/sqrt(disc);
  msg.Info()<<"     result (incl. flux): "<<m_res<<" +/- "<<m_error<<" ("<<m_error/m_res*100.<<" %)"<<endl;
} 


bool HD_PS_Base::WriteOut() {
  if ( m_read_dcfile ) {                // if DC file should be read
    system((string("mv \"")+m_path+m_file+string("\" \"")+m_path+m_file+string(".old\"")).c_str());

    ofstream to;
    to.open((m_path+m_file).c_str(),ios::out);

    // writes header
    to<<"# Decay: "<<p_hdc->ChannelName()<<endl;
    to<<"#        "<<p_hdc->ChannelNameNumbers()<<endl;
    // copy Channels, ME and Dalitz parameters
    char buffer[100];
    ifstream from;
    from.open((m_path+m_file+string(".old")).c_str());
    while (from.getline(buffer,100)) {
      if (buffer==string("Channels {")) {
        to<<"Channels {\n";
        from.getline(buffer,100);
        do {
          to<<buffer<<endl;
          from.getline(buffer,100);
        } while (buffer!=string("}")); 
        if (!m_foundPS) {                       // if there was no channel given
          for (int i=0;i<channels.size();i++) {
            if (channels[i]->Name()==string("Rambo") ||
                channels[i]->ChID()==string("Iso2"))
              to<<"    Isotropic"<<" "<<channels[i]->Alpha()<<";"<<endl;
            else
              to<<"    "<<channels[i]->ChID()<<" "<<channels[i]->Alpha()<<";"<<endl;
          }
        }
        to<<"}"<<endl;
      }
      if (buffer==string("Resonances {")) {
        to<<"Resonances {"<<endl;
        while (buffer!=string("}")) {
          from.getline(buffer,100);
          to<<buffer<<endl;
        }
      }
      if (buffer==string("Parameters {")) {
        to<<"Parameters {"<<endl;
        while (buffer!=string("}")) {
          from.getline(buffer,100);
          to<<buffer<<endl;
        }
      }
      if (buffer==string("Dalitz-Parameters {")) {
        to<<"Dalitz-Parameters {"<<endl;
        while (buffer!=string("}")) {
          from.getline(buffer,100);
          to<<buffer<<endl;
        }
      }
    }
    from.close();

    // write out result
    to<<"Result {"<<endl;
    to<<"   "<<m_res<<" "<<m_error<<" "<<m_max<<";"<<endl;
    to<<"}"<<endl;
    to.close();
  } // if (read DC file)
  else {                                // else create DC file                                  
    ofstream to;
    to.open((m_path+m_file).c_str(),ios::out);


    // writes header
    to<<"# Decay: "<<p_hdc->ChannelName()<<endl;
    to<<"#        "<<p_hdc->ChannelNameNumbers()<<endl;
    // write out channels
    to<<"Channels {\n"
      <<"\tAlwaysIntegrate = 0;   ! 0...read results and skip integration\n"
      <<"\t                       ! 1...don't read results and integrate\n";
    for (int i=0;i<channels.size();i++) {
      if (channels[i]->Name()==string("Rambo") ||
          channels[i]->ChID()==string("Iso2"))
        to<<"\tIsotropic"<<" "<<channels[i]->Alpha()<<";"<<endl;
      else
        to<<"\t"<<channels[i]->ChID()<<" "<<channels[i]->Alpha()<<";"<<endl;
    }
    to<<"}"<<endl;
    // write out result
    to<<"Result {"<<endl;
    to<<"   "<<m_res<<" "<<m_error<<" "<<m_max<<";"<<endl;
    to<<"}"<<endl;
    to.close();
  }
}

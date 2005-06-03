#include "HD_PS_Base.H"
#include "Hadron_Decay_Channel.H"
#include "Two_Body_PSs.H"
#include "Three_Body_PSs.H"
#include "Rambo.H"
#include "Data_Reader.H"
#include "Message.H"

using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

////////// class HD_Channel_Selector /////////

Single_Channel * HD_Channel_Selector::GetChannel( int nin, int nout, 
	const Flavour * flavs,string name )
{
  if (nin>1 || nout<2) {
    msg.Error()<<"Error in HD_Channel_Selector::GetChannel : "<<endl
	       <<"   No PS for channel ("<<nin<<" -> "<<nout<<" )"<<endl
	       <<"   Return nothing and hope for the best."<<endl;
    return NULL;
  }
  if (name==string("Isotropic")) return new Rambo(1,nout,flavs);
  if (nout==2) {
    if (name==string("Iso2")) return new Iso2Channel(flavs);
  }
  if (nout==3) {
    if (name==string("Dalitz_photon_23"))  return new Dalitz(flavs,Flavour(kf::photon),2,3);
    if (name==string("Dalitz_rho(770)_23")) return new Dalitz(flavs,Flavour(kf::rho_770),2,3);
  }

  msg.Error()<<"Error in HD_Channel_Selector::GetChannel : "<<endl
	     <<"   No channel for ("<<nin<<" -> "<<nout<<") with name "<<name<<endl
	     <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}

////////// class HD_PS_Base /////////

HD_PS_Base::HD_PS_Base( Hadron_Decay_Channel * hdc, vector<string> & _pst, 
	bool & mustinit, struct Model &_locmd ) :
  Multi_Channel(_pst[2]), p_hdc(hdc),
  p_channelselector(new HD_Channel_Selector), m_file(string("")),
  m_res(-1.), m_error(0.), m_max(-1.), m_flux(1./(2.*hdc->Flavours()[0].Mass()))
{
  if (_pst.size()>2) m_file = _pst[3]; 	// filename of DC file
  mustinit = Construct(_locmd);				// call Construct to do the rest
  delete p_channelselector;
}


HD_PS_Base::~HD_PS_Base() {}

void HD_PS_Base::Initialise() 
{
  CalculateNormalisedWidth();
  WriteOut();
}

bool HD_PS_Base::Construct(HADRONS::Model & _md )
{
  Initialize( _md );
  if (m_file!=string("")) {
    vector<vector<string> > helpsvv;
    Data_Reader reader = Data_Reader(string("|"),string(";"),string("!"));
    reader.SetAddCommandLine(false);
    reader.AddComment("#");
    reader.AddComment("//");
    reader.SetInputPath("./");
    reader.SetInputFile(m_file);
    reader.SetMatrixType(reader.MTransposed);
    if(!reader.MatrixFromFile(helpsvv)) {
      msg.Error()<<"ERROR in HD_PS_Base::Construct(HADRONS::Model&) :\n"
		 <<"   Read in failure, will abort."<<endl;
      abort();
    }

    string name;
    double weight;
    for (int i=0;i<helpsvv.size();i++) {
      if ( helpsvv[i][0]==string("Channels") ) {
		i++;											// next line
		while (helpsvv[i][0]!=string("}")) {
		  weight=1.;
		  if (helpsvv[i].size()>1) {					// if factor is given
			weight=atof(helpsvv[i][1].c_str());
		  }
		  AddChannel(helpsvv[i][0],weight);				// add decay channel
		  i++;
		}
      }
	  if ( helpsvv[i][0] == string("Parameters") ) {
		i++;
		while ( helpsvv[i][0] != string("}") ) {
		  cout<<i<<"  "<<helpsvv[i][0]<<"|"<<helpsvv[i][1]<<"|"<<helpsvv[i][2]<<endl;
		  if ( helpsvv[i][1] == string("=") ) {
			if ( helpsvv[i][0] == string("a") ) _md.pm.a   = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("b") ) _md.pm.b   = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("a2") ) _md.pm.a2 = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("b2") ) _md.pm.b2 = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("ME_MODEL") ) {
			  if ( string(helpsvv[i][2])=="Traces" )            _md.me = 1;
			  else if ( string(helpsvv[i][2])=="XYZ" )          _md.me = 2;
			  else if ( string(helpsvv[i][2])=="SimpleTraces" ) _md.me = 3;
			  else {
				msg.Error()<<"Error in HD_PS_Base::Construct ... "<<endl
				           <<"     There is no ME model with name "<<helpsvv[i][2]<<" implemented, yet."<<endl
						   <<"     Take XYZ instead."<<endl;
				_md.me = 2;
			  }
			}
			if ( string(helpsvv[i][0])=="RUNNING_WIDTH" ) _md.run = atoi( string(helpsvv[i][2]).c_str() );
			if ( string(helpsvv[i][0])=="FORM_FACTOR" )   _md.ff  = atoi( string(helpsvv[i][2]).c_str() );
		  }
		  i++;
		}
	  }
      if ( helpsvv[i][0]==string("Dalitz-Parameters") ) {
		i++;											// next line
		while (helpsvv[i][0]!=string("}")) {
		  if (helpsvv[i].size()==5) {					// there must be 5 given
			vector<double> dals;
			for (int j=0;j<5;j++) dals.push_back(atof(helpsvv[i][j].c_str()));
			//p_hdc->GetME()->SetDalitzParameters(dals);
		  }
		  i++;
		}
      }
      if (helpsvv[i][0]==string("Result")) {
		i++;
		while (helpsvv[i][0]!=string("}")) {
		  m_res   = atof(helpsvv[i][0].c_str());
		  m_error = atof(helpsvv[i][1].c_str());
		  m_max   = atof(helpsvv[i][2].c_str());
		  i++;
		}
		return false;
      }
    }
  }
  return true;
}

void HD_PS_Base::AddChannel(string name,double weight) {
  Single_Channel * sc = p_channelselector->GetChannel( 1, p_hdc->NOut(), 
	   											       p_hdc->Flavours(),
													   name );
  sc->SetAlpha(weight);
  Add(sc);									// add this to channels in Multi_Channel
}

void HD_PS_Base::CalculateNormalisedWidth() {
  msg_Tracking()<<"HD_PS_Base::CalculateNormalisedWidth()"<<endl;
  Reset();
  long int iter = Number()*5000*int(pow(2.,int(p_hdc->NOut())-2));
  int maxopt    = Number()*int(pow(2.,2*(int(p_hdc->NOut())-2)));
  
  long int n;
  int      opt=0;
  double   value, oldvalue=0., sum=0., sum2=0., result=-1., disc;
  bool     maxincrease, simple=false;
  while(opt<maxopt || m_error>0.01) {
    maxincrease = false;
    for (n=1;n<iter+1;n++) {
      value = p_hdc->Differential();
      sum  += value;
      sum2 += ATOOLS::sqr(value);
      AddPoint(value);
      if (value>m_max) { m_max = value; maxincrease = true; }
      if (value!=0. && value==oldvalue) { simple = true; break; }
      oldvalue = value;
    }
    opt++;
    Optimize(0.01);

    if (simple) break;
    n      = opt*iter;
    result = sum/n;
    disc   = sqr(sum/n)/((sum2/n - sqr(sum/n))/(n-1));
    if (disc>0) m_error  = result/disc;
  } 
  m_res  = m_flux*sum/n;
  m_max *= m_flux;
  disc   = sqr(m_res)/((sum2*sqr(m_flux)/n - sqr(m_res))/(n-1));
  if (disc>0) m_error  = result/disc;
}


bool HD_PS_Base::WriteOut() {
  system((string("mv ")+m_file+string(" ")+m_file+string(".old")).c_str());

  ofstream to;
  to.open(m_file.c_str(),ios::out);
  to<<"Channels {"<<endl;
  for (int i=0;i<channels.size();i++) {
    if (channels[i]->Name()==string("Rambo"))
      to<<"    Isotropic"<<" "<<channels[i]->Alpha()<<";"<<endl;
    else
      to<<"    "<<channels[i]->ChID()<<" "<<channels[i]->Alpha()<<";"<<endl;
  }
  to<<"}"<<endl;
  
  char buffer[100];
  ifstream from;
  from.open((m_file+string(".old")).c_str());
  while (from.getline(buffer,100)) {
    if (buffer==string("Dalitz-Parameters {")) {
      to<<"Dalitz-Parameters {"<<endl;
      while (buffer!=string("}")) {
	from.getline(buffer,100);
	to<<buffer<<endl;
      }
      break;
    }
  }
  from.close();
  to<<"Result {"<<endl;
  to<<"   "<<m_res<<" "<<m_error<<" "<<m_max<<";"<<endl;
  to<<"}"<<endl;
  to.close();
}

// Initialization of calculation model

void HD_PS_Base::Initialize( struct Model &md )
{
  md.me = 2;				// xyz functions
  md.run = 1;				// running width
  md.ff = 2;				// form factor
  md.pm.a   = 1.;
  md.pm.b	 = 1.;
  md.pm.a2  = 1.;
  md.pm.b2  = 1.;
  md.pm.GF  = 	1.16639e-5; 
  md.pm.Vud = 0.9744;
  md.pm.Vus = 0.2205;
  md.pm.fpi = 0.0924;
  md.pm.fK  = 0.113;
  md.pm.frho = 0.150;
  md.pm.grpp = 6.038;
}

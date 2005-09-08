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
	  ci.b=ab%10; ci.a=ab/10;	// int/int !
	  ci.nRes = 1;	
	}
	if( strcmp(result,"TwoResonances")==0 ) {
	  ci.name = result;
	  result = strtok( NULL, delim ); ci.res1 = result; 
	  result = strtok( NULL, delim ); ci.a = atoi( result );
	  result = strtok( NULL, delim ); ci.res2 = result; 
	  result = strtok( NULL, delim ); int bc = atoi( result );
	  ci.c=bc%10; ci.b=bc/10;	// int/int !
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
	  if( ci.res1==string("rho(770)") )  res.Set( kf::rho_770, md.pm.MR, md.pm.GR );
	  return new Dalitz(flavs,res,ci.a,ci.b);
	}
  }
  if (nout==4) {
	if( ci.name==string("TwoResonances") ) {
	  ResonanceFlavour res_a( ci.res1, md.pm.MA, md.pm.GA );
	  ResonanceFlavour res_v( ci.res2, md.pm.MV[0], md.pm.GR );
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
  Multi_Channel(_pst[2]), p_hdc(hdc),
  p_channelselector(new HD_Channel_Selector), m_file(string("")),
  m_res(-1.), m_error(0.), m_max(-1.), m_flux(1./(2.*hdc->Flavours()[0].Mass())),
  m_read_dcfile( read_dc ),
  m_path(_path), m_foundPS( false )
{
  if (_pst.size()>2) m_file = _pst[3]; 		// filename of DC file
  mustinit = Construct(_locmd);				// call Construct to do the rest
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
  Initialize( _md );
  if (m_file!=string("") && m_read_dcfile) {		// in case there is a DC file => read it !
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
		i++;											// next line
		while (helpsvv[i][0]!=string("}")) {
		  weight=1.;
		  if (helpsvv[i][0]==string("AlwaysIntegrate")) {
			skipresult = atoi( helpsvv[i][2].c_str() );
		  }
		  else {
			if (helpsvv[i].size()>1) {					// if factor is given
			  weight=atof(helpsvv[i][1].c_str());
			}
			else weight = 1.;
			if (IsChannel(helpsvv[i][0])) {				// if it is a channel that Sherpa can cope with
			  channels_vector.push_back( helpsvv[i][0] );		// save it for later
			  ch_weights.push_back( weight );
			}
		  }
		  i++;
		}
		m_foundPS = 1;								
		if (channels_vector.size()==0) {						// no channel found
		  AddChannel( string("Isotropic"), 1., _md );		//   take Rambo
		  skipresult = true;							//	 and don't read Result
		  m_foundPS = 0;				
		}
      }
	  if ( helpsvv[i][0] == string("Resonances") ) {
		i++;
		while (helpsvv[i][0]!=string("}")) {
		  if ( helpsvv[i][1] == string("->") ) {
			if ( helpsvv[i][0] == string("vector1") ) { _md.vec[0].i = atoi( string(helpsvv[i][2]).c_str() );
			                                            _md.vec[0].j = atoi( string(helpsvv[i][3]).c_str() ); }
			if ( helpsvv[i][0] == string("vector2") ) { _md.vec[1].i = atoi( string(helpsvv[i][2]).c_str() );
			                                            _md.vec[1].j = atoi( string(helpsvv[i][3]).c_str() ); }
			i++;
		  }
		}
	  }
	  if ( helpsvv[i][0] == string("Parameters") ) {
		i++;
		while ( helpsvv[i][0] != string("}") ) {
		  if ( helpsvv[i][1] == string("=") ) {
			// general constants
			if ( helpsvv[i][0] == string("V") ) _md.pm.a   = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("A") ) _md.pm.b   = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("V2") ) _md.pm.a2 = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("A2") ) _md.pm.b2 = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Vud")) _md.pm.Vud = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Vus")) _md.pm.Vus = atof( string(helpsvv[i][2]).c_str() );
			// weight factors for resonances
			if ( helpsvv[i][0] == string("alpha"))  _md.pm.alpha    = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("beta")) { _md.pm.beta     = atof( string(helpsvv[i][2]).c_str() );
													_md.pm.Beta[0]  = atof( string(helpsvv[i][2]).c_str() );
													_md.pm.Beta[1]  = atof( string(helpsvv[i][2]).c_str() ); }
			if ( helpsvv[i][0] == string("beta_1")) _md.pm.Beta[0]  = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("beta_2")) _md.pm.Beta[1]  = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("gamma"))  _md.pm.gamma    = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("delta"))  _md.pm.delta    = atof( string(helpsvv[i][2]).c_str() );
			// mass, width of resonances
			if ( helpsvv[i][0] == string("Mass_rho_770")) 	_md.pm.MR 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_rho_770")) 	_md.pm.GR 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_rho_1450")) 	_md.pm.MRR 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_rho_1450")) _md.pm.GRR 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_rho_1700")) 	_md.pm.MRRR = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_rho_1700")) _md.pm.GRRR = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_omega")) 	_md.pm.MO 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_omega")) 	_md.pm.GO 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_axial")) 		_md.pm.MA 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_axial"))		_md.pm.GA	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_vector"))  	{	_md.pm.MV[0]= atof( string(helpsvv[i][2]).c_str() );
			                                                	_md.pm.MV[1]= atof( string(helpsvv[i][2]).c_str() ); }
			if ( helpsvv[i][0] == string("Width_vector")) 	{	_md.pm.GV[0]= atof( string(helpsvv[i][2]).c_str() );
			                                                	_md.pm.GV[1]= atof( string(helpsvv[i][2]).c_str() ); }
			if ( helpsvv[i][0] == string("Mass_vector'")) 	{	_md.pm.Mv[0]= atof( string(helpsvv[i][2]).c_str() );
			                                                	_md.pm.Mv[1]= atof( string(helpsvv[i][2]).c_str() ); }
			if ( helpsvv[i][0] == string("Width_vector'"))	{	_md.pm.Gv[0]= atof( string(helpsvv[i][2]).c_str() );
			                                                	_md.pm.Gv[1]= atof( string(helpsvv[i][2]).c_str() ); }
			if ( helpsvv[i][0] == string("Mass_vector_1")) 	 	_md.pm.MV[0]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_vector_1"))	 	_md.pm.GV[0]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_vector_2")) 	 	_md.pm.MV[1]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_vector_2"))	 	_md.pm.GV[1]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_vector'_1"))  	_md.pm.Mv[0]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_vector'_1")) 	_md.pm.Gv[0]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Mass_vector'_2"))  	_md.pm.Mv[1]= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("Width_vector'_2")) 	_md.pm.Gv[1]= atof( string(helpsvv[i][2]).c_str() );
			// parameters for formfactor
			if ( helpsvv[i][0] == string("lambda")) _md.pm.lambda = atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("exp_alpha"))		_md.pm.exp_alpha 	= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("lambda0"))		_md.pm.lambda0		= atof( string(helpsvv[i][2]).c_str() );
			if ( helpsvv[i][0] == string("gamma_rho_770"))  _md.pm.gammaR		= atof( string(helpsvv[i][2]).c_str() );
			// further options
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
      if (helpsvv[i][0]==string("Result") && !skipresult ) {
		i++;
		while (helpsvv[i][0]!=string("}")) {
		  m_res   = atof(helpsvv[i][0].c_str());
		  m_error = atof(helpsvv[i][1].c_str());
		  m_max   = atof(helpsvv[i][2].c_str());
		  i++;
		}
		mustinit = false; 	// return: no need to integrate
      }
    }
	for (int i=0; i<channels_vector.size(); i++) {
	  AddChannel( channels_vector[i], ch_weights[i], _md ); 
	}
  }
  else {					// in case there is no DC file
	AddChannel( string("Isotropic"), 1., _md );
  }
  return mustinit;			// return: it has to be integrated or not
}

bool HD_PS_Base::AddChannel(string name,double weight,GeneralModel & md) 
{
  Single_Channel * sc = p_channelselector->GetChannel( 1, p_hdc->NOut(), 
	   											       p_hdc->Flavours(),
													   name, md );
  if (sc==NULL) return 0;
  sc->SetAlpha(weight);
  Add(sc);									// add this to channels in Multi_Channel
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
  while(opt<maxopt || m_error>0.01) {
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

    if (simple) break;
    n      = opt*iter;
    result = sum/n;
    disc   = sqr(sum/n)/((sum2/n - sqr(sum/n))/(n-1));
    if (disc>0) m_error  = result/sqrt(disc);
	msg.Info()<<"     result (w/o flux): "<<result<<" +/-"<<m_error<<" ("<<m_error/result*100.<<" %)"<<endl;
  } 
  m_res  = m_flux*sum/n;
  m_max *= m_flux;
  disc   = sqr(m_res)/((sum2*sqr(m_flux)/n - sqr(m_res))/(n-1));
  if (disc>0) m_error  = result/sqrt(disc);
} 


bool HD_PS_Base::WriteOut() {
  if ( m_read_dcfile ) {				// if DC file should be read
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
		if (!m_foundPS) { 						// if there was no channel given
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
  else {								// else create DC file									
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

// Initialization of calculation model

void HD_PS_Base::Initialize( struct GeneralModel &md )
{
  md.me = 2;				// xyz functions
  md.run = 3;				// running width
  md.ff = 1;				// form factor
  md.vec[0].i = 1;			// number of decay product of vector resonance 1
  md.vec[0].j = 3;			// number of decay product of vector resonance 1
  md.vec[1].i = 2;			// number of decay product of vector resonance 2
  md.vec[1].j = 3;			// number of decay product of vector resonance 2
  md.pm.a   = 1.;
  md.pm.b	= 1.;
  md.pm.a2  = 1.;
  md.pm.b2  = 1.;
  md.pm.GF  = 1.16639e-5; 
  md.pm.Vud = 0.9744;
  md.pm.Vus = 0.2205;
  md.pm.fpi = 0.0924;
  md.pm.fK  = 0.113;
  md.pm.frho = 0.150;
  md.pm.grpp = 6.038;
  md.pm.MR   = 0.775;
  md.pm.GR   = 0.14;
  md.pm.MO	 = 0.78;
  md.pm.GO   = 0.005;
  md.pm.MRR  = 1.465;
  md.pm.GRR  = 0.4;
  md.pm.MRRR = 1.72;
  md.pm.GRRR = 0.2;
  md.pm.alpha = 0.;
  md.pm.beta  = 0.;
  md.pm.gamma = 0.;
  md.pm.delta = 0.;
  md.pm.Beta[0] = md.pm.beta;
  md.pm.Beta[1] = md.pm.beta;
  md.pm.MV[0] = md.pm.MR;
  md.pm.MV[1] = md.pm.MR;
  md.pm.GV[0] = md.pm.GR;
  md.pm.GV[1] = md.pm.GR;
  md.pm.Mv[0] = md.pm.MRR;
  md.pm.Mv[1] = md.pm.MRR;
  md.pm.Gv[0] = md.pm.GRR;
  md.pm.Gv[1] = md.pm.GRR;
  md.pm.MA    = 1.205;
  md.pm.GA    = 0.44;
}

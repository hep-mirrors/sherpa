#include "Decay_Table_Reader.H"
#include "Hadron_Decay_Channel.H"
#include "Decay_Table.H"
#include "Data_Reader.H"
#include "Message.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Decay_Table_Reader::Decay_Table_Reader(string path,string file) :
  m_path( path ), m_file( file ), m_fulldecay(1), m_antidecays(1), m_createbooklet(0)
{
  if (file==string("")) {
    msg.Error()<<"Error in Decay_Table_Reader::Decay_Table_Reader("
	       <<path<<","<<file<<") : "<<endl
	       <<"   No file specified, will return and hope for the best."<<endl;
    return;
  }
  msg_Info()<<"New Decay_Table_Reader("<<path+file<<")"<<endl;
  Data_Reader reader = Data_Reader(string("XxX"),string(";"),string("!")); // only split b/w words
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  reader.SetMatrixType(mtc::transposed);
  // read decay table file
  vector<vector<string> > helpsvv;
  if(!reader.MatrixFromFile(helpsvv)) {
    msg.Error()<<"ERROR in Decay_Table_Reader::"
	       <<"Decay_Table_Reader(string,string) :\n"
	       <<"   Read in failure "<<path<<file<<", will abort."<<endl;
    abort();
  }
  // split b/w |'s : three col's 
  m_helpsvv.clear();
  vector<string> line;
  bool finishedBR (false);
  for( size_t row=0; row<helpsvv.size(); ++row ) {
    line.clear();
    int col=0;
    // products, BR info, DC filename (optional)
    for( int info_col=0; info_col<3; ++info_col ) {
      string entry ("");
      if( col<helpsvv[row].size() ) entry = helpsvv[row][col];
      col++;
      while( col<helpsvv[row].size() && helpsvv[row][col]!=string("|") ) {
        entry += string(" ")+helpsvv[row][col];
        col++;
      }
      if( entry != string("") ) line.push_back( entry );
      if( col==helpsvv[row].size() ) break;
      col++;
    }
    m_helpsvv.push_back(line);
  }
}

int Decay_Table_Reader::FillDecayTable(Decay_Table * dt)
{
  int              nchannels = 0;									// number of channels
  vector<int>      helpkfc;											// aux. kf-code used to extract flavours
  double           BR, dBR;											// branching ratio with error
  string           sbr, sdbr;                                       // string version of BR and dBR
  string           origin;                                          // "origin" of channel (PDG,intuition,...)
  Flavour          flav;											// involved flavors
  Decay_Channel  * dc;												// pointer on decay channel

  for (int i=0;i<m_helpsvv.size();i++) {							// for each element
	if (m_helpsvv[i][0] != string("NO_FULLDECAY") &&				// if no special keyword
	    m_helpsvv[i][0] != string("NO_ANTI") &&
		m_helpsvv[i][0] != string("CREATE_BOOKLET") ) {
	  if (ExtractFlavours(helpkfc,m_helpsvv[i][0])) {				//   extract flavours
        ExtractBRInfo( m_helpsvv[i][1], sbr, sdbr, origin );
		BR = double(atof(sbr.c_str()));					            //   get branching ratio
		dBR = double(atof(sdbr.c_str()));				            //   get branching ratio error
        if ( BR<1.e-6 ) BR=1.e-6;
        dc = new Decay_Channel(dt->Flav());                         //   create new decay channel
        for (int j=0;j<helpkfc.size();++j) {
          flav = Flavour(kf::code(abs(helpkfc[j])));
          if (helpkfc[j]<0) flav = flav.Bar();
          dc->AddDecayProduct(flav);				
        }
        dc->SetWidth(BR*dt->Flav().Width());						//   set decay width
        dc->SetDeltaWidth(dBR*dt->Flav().Width());					//   set decay width
        dc->SetOrigin(origin);                                      //   set origin
        dc->SetProcessName();										//   set process name
        dc->SetMEType("");
        if (m_helpsvv[i].size()>2) dc->SetPSFile(m_helpsvv[i][2]);
        dt->AddDecayChannel(dc);									//   add to decay channels
        nchannels++;								
	  }
	}
	else {															// special keywords
	  if (m_helpsvv[i][0] == string("NO_FULLDECAY"))	m_fulldecay  = 0;
	  if (m_helpsvv[i][0] == string("NO_ANTI")) 		m_antidecays = 0;
      m_createbooklet = 0;
    }
  } // for
  return nchannels;												// return number of channels
}

void Decay_Table_Reader::FillInMatrixElementsAndPS( Decay_Table * dt, Channel_Map * chmap, GeneralModel startmd )
{
  Decay_Channel        * dc;								// pointer on decay channel
  Hadron_Decay_Channel * hdc;								// pointer on hadronic decay channel
  int rewrite (0);										// rewrite H file ?
  for (int i=0;i<dt->NumberOfDecayChannels();i++) {			// for each channel in decay table
    dc  = dt->GetDecayChannel(i);							//   get channel
    hdc = new Hadron_Decay_Channel(dc,m_path);				//   create new hadronic channel
    rewrite += hdc->InitialisePhaseSpace(m_helpsvv[i],startmd);		//   initialise phase space
	hdc->SetFullDecay((m_fulldecay<<1)+m_antidecays);		//   fulldecay? antidecay?
	hdc->SetCreateBooklet( m_createbooklet );				// 	 create booklet?
	(*chmap)[ dc ] = hdc;									//   map channel to hadronic channel
  }
  if (rewrite) {											// rewrite H file
	system( (string("mv \"")+m_path+m_file+string("\" \"")+m_path+m_file+string(".old\"")).c_str() );
	ofstream f( (m_path + m_file).c_str() );
	f<<"# outgoing part. \t | BR(Delta BR) \t [Origin] \t | DC-file"<<endl;
	for (int j=0; j<dt->NumberOfDecayChannels();j++) {
	  int nr (m_helpsvv[j].size());							// # columns
	  for (int i=0; i<nr; i++) {
		f<<m_helpsvv[j][i];
		if (i==nr-1) f<<";"<<endl;
		else f<<"\t | ";
	  }
	}
	f<<endl;
	if (hdc->FullDecay() & 2) f<<"! ";
	f<<"NO_FULLDECAY   ! use this option: no decay of daughters\n";
	if (hdc->FullDecay() & 1) f<<"! ";
	f<<"NO_ANTI        ! use this option: anti-mother does not decay\n";
	f.close();
  }
}

bool Decay_Table_Reader::ExtractFlavours(vector<int> & helpkfc,string help)
{
  helpkfc.clear();													
  size_t pos = help.find("{");
  bool             hit;
  if (pos!=string::npos) help = help.substr(pos+1);
  else {
    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
	       <<"   Something wrong with final state of decay (Bracket missing) :"<<help<<endl
           <<"   Please check "<<m_path<<m_file<<"!"<<endl
	       <<"   Will skip it."<<endl;
    return false;
  }
  pos    = help.find("}");
  if (pos!=string::npos) help = help.substr(0,pos);
  else {
    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
	       <<"   Something wrong with final state of decay (Bracket missing) :"<<help<<endl
           <<"   Please check "<<m_path<<m_file<<"!"<<endl
	       <<"   Will skip it."<<endl;
    return false;
  }
  hit    = true;
  while (hit) {
    pos      = help.find(",");
    if (pos!=string::npos) {
      helpkfc.push_back(atoi((help.substr(0,pos)).c_str()));
      help  = help.substr(pos+1);
    }
    else {
      helpkfc.push_back(atoi(help.c_str()));
      hit = false;
    }
  }
//  if (helpkfc.size()<2) {
//    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
//	       <<"   Something wrong with final state of decay (Too little particles) : ";
//    for (int j=0;j<helpkfc.size();j++) msg.Error()<<helpkfc[j]<<" ";
//    msg.Error()<<endl<<"   Will skip it."<<endl;
//    return false;
//  } 
  if (helpkfc.size()<1) {
    msg.Error()<<"WARNING in Decay_Table_Reader:: : "<<endl
	       <<"   Something wrong with final state of decay. (no particles?)"<<endl
           <<"   Please check "<<m_path<<m_file<<"!"<<endl 
           <<"   Will skip it and hope for the best."<<endl;
    return false;
  } 
  return true;
}
 
void Decay_Table_Reader::ExtractBRInfo( string entry, string & sbr, string & sdbr, string & sorigin )
{
  size_t posa, posb;        // start and end of things b/w brackets
  size_t posmin;            // start of first bracket
  // extract Delta BR
  posa = entry.find("(");
  posb = entry.find(")");
  posmin = posa;
  if( posa!=string::npos && posb!=string::npos ) sdbr = entry.substr(posa+1,posb-posa-1);
  else                                           sdbr = string("-1");         // no Delta BR given
  // extract Origin
  posa = entry.find("[");
  posb = entry.find("]");
  if( posmin==string::npos || (posmin!=string::npos && posmin>posa) ) posmin = posa;
  if( posa!=string::npos && posb!=string::npos ) sorigin = entry.substr(posa+1,posb-posa-1);
  else                                           sorigin = string("PDG");   // no Origin
  // extract BR
  if( posmin!=string::npos ) sbr = entry.substr(0,posmin);
  else                       sbr = entry.substr(0);
}


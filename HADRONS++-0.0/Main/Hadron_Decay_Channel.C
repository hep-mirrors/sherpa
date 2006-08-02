#include "Hadron_Decay_Channel.H"
#include "HD_ME_Base.H"
#include "Spin_Correlation_Tensor.H"
#include "Poincare.H"
#include "Random.H"
#include "XYZFuncs.H"
#include "MyStrStream.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel( Decay_Channel * _dc, string _path ) :
    p_dc(_dc),
  Integrable_Base(1,_dc->NumberOfDecayProducts()),
  m_path(_path), m_fulldecay(3), m_mass_smearing(true), m_createbooklet(0),
  p_indices(NULL), p_ampls(NULL)
{
  m_resultpath      = string("./");                             // where to write results
  m_resultfile      = m_histofile = string("");                 // filename
  m_name            = p_dc->ProcessName();                      // name of proces
  p_flavours        = new Flavour[m_nin+m_nout];                    
  p_flavours[0]     = p_dc->GetDecaying();                      // decaying particle
  m_channelname     = string("");
  m_chnamenumbers   = string("0");
  m_channelname     = p_flavours[0].IDName() + string(" --> ");
  m_chnamenumbers.append( (p_flavours[0].IDName()).length()-1, ' ' );
  m_chnamenumbers.append(" --> ");
  char helpch[2];
  for (int i=0;i<m_nout;i++) {           // decay products
    p_flavours[i+1] = p_dc->GetDecayProduct(i);
    m_channelname  += p_flavours[i+1].IDName() + string(" ");
    sprintf( helpch, "%i%", i+1 );
    m_chnamenumbers.append(string(helpch));
    m_chnamenumbers.append( (p_flavours[i+1].IDName()).length(), ' ' );
  }
  
  p_indices = new vector<pair<int,int> >;                        // index bookkeeping
  p_ampls   = new vector<Complex>;                              // new amplitude tensor

  p_ps = new HD_PS_Base(this/*,m_path,PStype,mustinit,locmd,selectionlines,read_dc*/);

//   HD_ME_Selector mesel;                                         // ME selector
//   p_me = mesel.GetME(m_nin,m_nout,p_flavours);                  // get the appropr. ME
//   p_me->SetPath( m_path );                                      // set Decaydata path
//   msg.Tracking()<<"Matrix Element for "<<m_channelname<<" : "<<p_me->METype()<<"."<<endl;
//   msg.Tracking()<<"Currents for "<<m_channelname<<endl;     // new for currents
//   p_currents = mesel.GetCurrents(m_nin,m_nout,p_flavours,m_mefactor);      // new for currents
//   if(p_currents[0]) p_currents[0]->SetPath( m_path );
//   if(p_currents[1]) p_currents[1]->SetPath( m_path );
  // check for identical particles
  Flavour refflav;
  double symfactor (1);         
  int l(0), lfac (1);              
  for( int i=0; i<m_nout; ++i ) {
    refflav = p_flavours[i+1];
    l = 0;
    lfac = 1;
    for( int j=0; j<m_nout; ++j ) {
      if( p_flavours[j+1]==refflav ) {
        l++;
        lfac *= l;
      }
    }
    symfactor *= lfac;
  }
  m_symmetry = 1./sqrt(symfactor);
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
  if (p_dc) { delete p_dc; p_dc=NULL; }
  if (p_ps) { delete p_ps; p_ps=NULL; }
//   if (p_me) { delete p_me; p_me=NULL; }
//   if (p_currents) { // new for currents
//     if(p_currents[0]) delete p_currents[0]; p_currents[0]=NULL; 
//     if(p_currents[1]) delete p_currents[1]; p_currents[1]=NULL; 
//     delete[] p_currents; p_currents=NULL;
//   }
  if (p_ampls) { delete p_ampls; p_ampls=NULL; }
  if (p_indices) { delete p_indices; p_indices=NULL; }
}


bool Hadron_Decay_Channel::Initialise(vector<string> & PStype, GeneralModel startmd)
  // PStype: decay products | BR | DC filename  <-- line of Decay file
{
  m_startmd=startmd;
  bool mustinit;
  bool rewriteDT (false);
  if ( PStype.size() == 2 ) {                                   // in case no DC file given
    string fn("");                                              // filename of DC file
    fn += p_dc->GetDecaying().ShellName() + string("_");
    for ( int i=0; i<p_dc->NumberOfDecayProducts(); i++ ) {
      fn += p_dc->GetDecayProduct(i).ShellName();
    }
    fn += string(".dat");
    PStype.push_back( fn );                                     // generate DC filename
    rewriteDT = true;                                            // rewrite hadron decay file
  }
  // check if dc file exists
  m_filename = PStype[2];
  ifstream dcf( (m_path+m_filename).c_str() );
  bool read_dc = dcf;                               // read DC file if it exists
  dcf.close();

  if (read_dc) {
    msg.Tracking()<<METHOD<<": read "<<m_path<<m_filename<<endl;
    Data_Reader reader(string("|"),string(";"),string("!"));
    reader.SetAddCommandLine(false);
    reader.AddComment("#");
    reader.AddComment("//");
    reader.SetInputPath(m_path);
    reader.SetInputFile(m_filename);
    reader.SetMatrixType(mtc::transposed);

    // process <Options>
    vector<vector<string> > options_svv;
    reader.SetFileBegin("<Options>"); reader.SetFileEnd("</Options>");
    if(reader.MatrixFromFile(options_svv)) ProcessOptions(options_svv);
    else {
      msg.Error()<<METHOD<<": Error.\n"
         <<"  Read in failure for <Options> section in "<<m_path<<m_filename<<", will abort."<<endl;
      abort();
    }

    // process <ME>
    vector<vector<string> > me_svv;
    GeneralModel model_for_ps;
    reader.SetFileBegin("<ME>"); reader.SetFileEnd("</ME>");
    reader.RereadInFile();
    if(reader.MatrixFromFile(me_svv)) ProcessME(me_svv, reader, model_for_ps);
    else {
      msg.Error()<<METHOD<<": Error.\n"
                 <<"  Read in failure for <Phasespace> section in "<<m_path
                 <<m_filename<<", will abort."<<endl;
      abort();
    }

    // process <Phasespace>
    vector<vector<string> > ps_svv;
    reader.SetFileBegin("<Phasespace>"); reader.SetFileEnd("</Phasespace>");
    reader.RereadInFile();
    if(!reader.MatrixFromFile(ps_svv)) {
    msg.Error()<<METHOD<<": Error."<<"  Read in failure for <Phasespace> section in "
        <<m_path<<m_filename<<", will abort."<<endl;
      abort();
    }
    ProcessPhasespace(ps_svv, reader, model_for_ps);

    // process <Result> // don't do it before ME and phasespace, or CalcNormWidth doesn't work!
    vector<vector<string> > result_svv;
    reader.SetFileBegin("<Result>"); reader.SetFileEnd("</Result>");
    reader.RereadInFile();
    reader.MatrixFromFile(result_svv);
    ProcessResult(result_svv);
  }
  else { // if DC file does not exist yet
    m_mes.push_back( MEPair(1.0,new Isotropic(m_nout,p_flavours)) );
    p_ps->AddChannel( string("Isotropic"), 1., m_startmd );
    vector<double> results = p_ps->CalculateNormalisedWidth();
    SetAlwaysIntegrate(false);
    SetMassSmearing(false);
    WriteOut(results, true);
  }
  return rewriteDT;
}

void Hadron_Decay_Channel::ProcessOptions(vector<vector<string> > helpsvv)
{
  for (int i=0;i<helpsvv.size();i++) {
    if (helpsvv[i][0]==string("AlwaysIntegrate")) {
      SetAlwaysIntegrate( atoi( helpsvv[i][2].c_str() ) );
    }
    else if(helpsvv[i][0]==string("MassSmearing")) {
      SetMassSmearing( atoi(helpsvv[i][2].c_str()) );
    }
  }
}

void Hadron_Decay_Channel::ProcessPhasespace(vector<vector<string> > ps_svv,
                                             Data_Reader           & reader,
                                             GeneralModel const    & model_for_ps)
{
  int nr_of_channels=0;
  for (int i=0;i<ps_svv.size();i++) {
    double weight = ToType<double>(ps_svv[i][0]);
    if( p_ps->AddChannel( ps_svv[i][1], weight, model_for_ps ) ) nr_of_channels++;
    else {
      msg.Error()<<METHOD<<": Warning. "<<ps_svv[i][1]<<"in"<<m_path<<m_filename
        <<" is not a valid phase space channel. Will ignore it."<<endl;
    }
  }
  if(nr_of_channels == 0) {
    msg.Error()<<METHOD<<": Warning. No valid phase space channels found in "
      <<m_path<<m_filename<<". Using Isotropic."<<endl;
    p_ps->AddChannel( string("Isotropic"), 1., m_startmd );
  }
}

void Hadron_Decay_Channel::ProcessME( vector<vector<string> > me_svv,
                                      Data_Reader           & reader,
                                      GeneralModel          & model_for_ps )
{
  int nr_of_mes=0;
  for (int i=0;i<me_svv.size();i++) {
    double factor = ToType<double>(me_svv[i][0]);
    if(me_svv[i].size()==2) {
      msg.Tracking()<<"Selecting ME for "<<m_channelname<<endl;
      HD_ME_Base* me = SelectME( me_svv[i][1] );
      me->SetPath(m_path);
      msg.Tracking()<<"  "<<me->METype()<<endl;
      vector<vector<string> > parameter_svv;
      reader.SetFileBegin("<"+me_svv[i][1]+">"); reader.SetFileEnd("</"+me_svv[i][1]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter_svv);
      GeneralModel me_model = Parameters2Model(parameter_svv, &model_for_ps);
      me->SetModelParameters( me_model );
      m_mes.push_back( MEPair(ToType<double>(me_svv[i][0]),me) );
      nr_of_mes++;
    }
    if(me_svv[i].size()==3) {
      msg.Tracking()<<"Selecting currents for "<<m_channelname<<endl;
      Current_Base* current1 = SelectCurrent(me_svv[i][1]);
      current1->SetPath(m_path);
      vector<vector<string> > parameter1_svv;
      reader.SetFileBegin("<"+me_svv[i][1]+">"); reader.SetFileEnd("</"+me_svv[i][1]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter1_svv);
      GeneralModel current1_model = Parameters2Model(parameter1_svv, &model_for_ps);
      current1->SetModelParameters( current1_model );

      Current_Base* current2 = SelectCurrent(me_svv[i][2]);
      current2->SetPath(m_path);
      vector<vector<string> > parameter2_svv;
      reader.SetFileBegin("<"+me_svv[i][2]+">"); reader.SetFileEnd("</"+me_svv[i][2]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter2_svv);
      GeneralModel current2_model = Parameters2Model(parameter2_svv, &model_for_ps);
      current2->SetModelParameters( current2_model );

      msg.Tracking()<<"  "<<current1->m_name<<endl;
      msg.Tracking()<<"  "<<current2->m_name<<endl;

      // Sanity checks for current selection
      if( 1+m_nout != current1->m_n + current2->m_n ) {
        msg.Error()<<"Error in "<<METHOD<<": Current selection does not look sane "
            <<"for "<<m_channelname<<". Check decaychannelfile."<<std::endl;
        abort();
      }

      m_currents.push_back( CurrentsPair( ToType<double>(me_svv[i][0]),
                                          pair<Current_Base*,Current_Base*>(current1,current2) ) );
      nr_of_mes++;
    }
  }
  if(nr_of_mes == 0) {
  msg.Error()<<METHOD<<": Warning. No valid matrix element found in "
      <<m_path<<m_filename<<". Using Isotropic."<<endl;
    m_mes.push_back( MEPair(1.0,new Isotropic(m_nout,p_flavours)) );
  }
}

void Hadron_Decay_Channel::ProcessResult(vector<vector<string> > result_svv)
{
  if(result_svv.size()!=1 || m_always_integrate) {
    vector<double> results = p_ps->CalculateNormalisedWidth();
    WriteOut(results);
  }
  else {
    p_ps->SetResult(ToType<double>(result_svv[0][0]));
    p_ps->SetError(ToType<double>(result_svv[0][1]));
    p_ps->SetMaximum(ToType<double>(result_svv[0][2]));
  }
}

GeneralModel Hadron_Decay_Channel::Parameters2Model(vector<vector<string> > helpsvv,
                                                    GeneralModel          * other_model )
{
  GeneralModel model(m_startmd);
  Algebra_Interpreter ip;
  for (int i=0;i<helpsvv.size();i++) {
    if ( helpsvv[i][1] == string("=")) {
      if( helpsvv[i].size() == 3 ) {        // <name> = <real value>
        double real = ToType<double> (ip.Interprete(helpsvv[i][2]) );
        model[helpsvv[i][0]] = real;
        if(other_model) {
          (*other_model)[helpsvv[i][0]] = real;
        }
      }
      if( helpsvv[i].size() == 4 ) {        // <name> = <complex value>
        double abs   = ToType<double>(ip.Interprete(helpsvv[i][2]) );
        double phase = ToType<double>( ip.Interprete(helpsvv[i][3]) );
        model[helpsvv[i][0]+string("_abs")] = abs;
        model[helpsvv[i][0]+string("_phase")] = phase;
        if(other_model) {
          (*other_model)[helpsvv[i][0]+string("_abs")] = abs;
          (*other_model)[helpsvv[i][0]+string("_phase")] = phase;
        }
      }
    }
    if ( helpsvv[i][2] == string("=")) {
      if( helpsvv[i].size() == 4 ) {        // <name> <index> = <real value>
        double real = ToType<double>( ip.Interprete(helpsvv[i][3]) );
        model[helpsvv[i][0]+string("_")+helpsvv[i][1]] = real;
        if(other_model) (*other_model)[helpsvv[i][0]+string("_")+helpsvv[i][1]] = real;
      }
      if( helpsvv[i].size() == 5 ) {        // <name> <index> = <complex value>
        double abs   = ToType<double> ( ip.Interprete(helpsvv[i][3]) );
        double phase = ToType<double> ( ip.Interprete(helpsvv[i][4]) );
        model[helpsvv[i][0]+string("_")+helpsvv[i][1]+string("_abs")] = abs;
        model[helpsvv[i][0]+string("_")+helpsvv[i][1]+string("_phase")] = phase;
        if(other_model) {
          (*other_model)[helpsvv[i][0]+string("_")+helpsvv[i][1]+string("_abs")] = abs;
          (*other_model)[helpsvv[i][0]+string("_")+helpsvv[i][1]+string("_phase")] = phase;
        }
      }
    }
  }
  return model;
}

void Hadron_Decay_Channel::WriteModelOnScreen( GeneralModel _locmd )
{
  if( msg.LevelIsDebugging() ) {
  msg.Out()
    <<"-----------------------------------------------------\n"
    <<"Modelparameters for channel: "<<m_channelname<<"\n"<<endl;
  GeneralModel::iterator md_it;
  for ( md_it = _locmd.begin(); md_it != _locmd.end(); ++md_it ) {
    msg.Out()<<"   "<<md_it->first<<":\t"<<md_it->second<<endl;
  }
  msg.Out()
    <<"-----------------------------------------------------"<<endl;
  }
}

bool Hadron_Decay_Channel::WriteOut( vector<double> results, bool newfile ) {

  if ( newfile ) {                // if DC file doesn't exist yet
    ofstream to;
    to.open((m_path+m_filename).c_str(),ios::out);

    // write header
    to<<"# Decay: "<<ChannelName()<<endl;
    to<<"#        "<<ChannelNameNumbers()<<endl<<endl;

    // write out options
    to<<"<Options>"<<endl;
    to<<"  AlwaysIntegrate = "<<m_always_integrate<<"    # 0...read results and skip integration"<<endl;
    to<<"                         # 1...don't read results and integrate"<<endl;
    to<<"  MassSmearing    = "<<m_mass_smearing<<"    # 0...turn off mass smearing for this decay's products"<<endl;
    to<<"                         # 1...turn on mass smearing for this decay's products"<<endl;
    to<<"</Options>"<<endl<<endl;

    // write out phasespace settings
    to<<"<Phasespace>"<<endl;
    to<<"  1.0 Isotropic"<<endl;
    to<<"</Phasespace>"<<endl<<endl;

    // write out ME settings
    to<<"<ME>"<<endl;
    to<<"  1.0 Isotropic"<<endl;
    to<<"</ME>"<<endl<<endl;

    // write out result
    to<<"<Result>"<<endl;
    to<<"  "<<results[0]<<" "<<results[1]<<" "<<results[2]<<";"<<endl;
    to<<"</Result>"<<endl;
    to.close();
  } // if (read DC file)
  else {                                // if DC file exists
    system((string("mv \"")+m_path+m_filename+string("\" \"")+m_path+"."+m_filename+string(".old\"")).c_str());
    ofstream to((m_path+m_filename).c_str(),ios::out);

    // copy Options, Phasespace, ME, ...
    char buffer[100];
    ifstream from;
    from.open((m_path+"."+m_filename+string(".old")).c_str());
    while (from.getline(buffer,100)) {
      if (buffer==string("<Result>")) break;
      else to<<buffer<<endl;
    }
    from.close();

    // write out result
    to<<endl;
    to<<"<Result>"<<endl;
    to<<"  "<<results[0]<<" "<<results[1]<<" "<<results[2]<<";"<<endl;
    to<<"</Result>"<<endl;
    to.close();
  }
  return 1;
}


// w/o spin correlation
double Hadron_Decay_Channel::Differential()
{
  p_momenta[0] = Vec4D(p_flavours[0].PSMass(),0.,0.,0.);   // decay from rest
  p_ps->GeneratePoint(p_momenta);                               // generate a PS point
  p_ps->GenerateWeight(p_momenta);                              // calculate its weight factor
  double weight = p_ps->Weight();                               // get weight factor

  vector<Complex>* ampls = new vector<Complex>;
  vector<pair<int,int> >* indices = new vector<pair<int,int> >;

  if(m_mes.size()>0) {
    // get the amplitude vector for the first ME
    m_mes[0].second->operator()(  p_momenta, p_ampls, p_indices, 1 );
    // add the amplitudes for each other ME
    vector<MEPair>::iterator mpit;
    vector<MEPair>::iterator mpstart = m_mes.begin();
    mpstart++;
    for(mpit=mpstart; mpit!=m_mes.end();mpit++) {
      ((*mpit).second)->operator()(  p_momenta, ampls, indices, 1 );
      AddAmpls(ampls,indices);
      ampls->clear(); indices->clear();
    }
  }

  if(m_currents.size()>0) {
    // get the amplitude vector for the first currents (only if not initialised yet)
    if( m_mes.size()==0 )
      ContractCurrents(m_currents[0].second.first, m_currents[0].second.second,
                       p_momenta, 1, m_currents[0].first,
                       p_ampls, p_indices);
    // add the amplitudes for all (further) currents
    vector<CurrentsPair>::iterator cpit;
    vector<CurrentsPair>::iterator cpstart = m_currents.begin();
    if( m_mes.size()==0 ) cpstart++;
    for(cpit=cpstart; cpit!=m_currents.end();cpit++) {
      ContractCurrents((*cpit).second.first, (*cpit).second.second,
                      p_momenta, 1, (*cpit).first,
                      ampls, indices);
      AddAmpls(ampls,indices);
      ampls->clear(); indices->clear();
    }
  }

  delete ampls; ampls=NULL; delete indices; indices=NULL;

//   if(p_currents && p_currents[0] && p_currents[1]) {
//     ContractCurrents(p_currents[0], p_currents[1],
//                      p_momenta,1,m_mefactor,
//                      p_ampls,p_indices);
//   }
//   else {
//     (*p_me)(  p_momenta,                                          // phase space point
//               p_ampls, p_indices,                                 // ampl. tensor and indices
//               1 );                                            // spinor base
//   }
  double value (0.);
  if( p_ampls->size() ) {
    for( size_t i=0; i<p_ampls->size(); ++i ) {
      value += norm( (*p_ampls)[i] );
    }
    value /= (p_dc->GetDecaying().IntSpin()+1.);
    if(!(value>=0.0 || value<0.0)) {
      msg.Error()<<"Problem in Differential: ME or currents delivered a nan:"<<std::endl;
//       if(p_currents[0]) msg.Error()<<"  current[0]: "<<*(p_currents[0])<<std::endl;
//       if(p_currents[1]) msg.Error()<<"  current[1]: "<<*(p_currents[1])<<std::endl;
//       msg.Error()<<"  Amplitudes:"<<std::endl;
//       for(int i=0; i<p_ampls->size(); i++) {
//         msg.Error()<<"    "<<i<<": "<<(*p_ampls)[i]<<std::endl;
//       }
      msg.Error()<<"value="<<value<<std::endl;
      abort();
    }
  }
  else value = 1.;                                              // isotropic
  return value*weight*m_symmetry;
}


// with spin correlation
double Hadron_Decay_Channel::Differential( Vec4D * mom, Spin_Density_Matrix * sigma )
{
  if( !mom ) return Differential();                             // if no momentum
  double ret;
  if( !sigma ) {                                                // no SDM
    ret = Differential();
    // boost into Lab frame
    Poincare lambda(mom[0]);
    lambda.Invert();
    for( int i=0; i<m_nout+1; ++i ) {
      p_momenta[i] = lambda*p_momenta[i];
      mom[i] = p_momenta[i];
    }
    return ret;
  }
  // spin correlations
  // get PS point in rest frame
  p_momenta[0] = Vec4D(p_flavours[0].PSMass(),0.,0.,0.);        // decay from rest
  p_ps->GeneratePoint(p_momenta);                               // generate a PS point
  p_ps->GenerateWeight(p_momenta);                              // calculate its weight factor
  double weight = p_ps->Weight();                               // weight factor
  // boost into Lab system
  Poincare lambda(mom[0]);
  lambda.Invert();
  for( int i=0; i<m_nout+1; ++i ) {
    p_momenta[i] = lambda*p_momenta[i];
    mom[i] = p_momenta[i];
  }

  // get amplitude tensor
  vector<Complex>* ampls = new vector<Complex>;
  vector<pair<int,int> >* indices = new vector<pair<int,int> >;
  
  vector<MEPair>::iterator mpit;
  for(mpit=m_mes.begin(); mpit!=m_mes.end();mpit++) {
    (*mpit).second->operator()(  p_momenta, ampls, indices, 1 );
    AddAmpls(ampls,indices);
    ampls->clear(); indices->clear();
  }
  vector<CurrentsPair>::iterator cpit;
  for(cpit=m_currents.begin(); cpit!=m_currents.end();cpit++) {
    ContractCurrents((*cpit).second.first, (*cpit).second.second,
                     p_momenta, 1, (*cpit).first,
                     ampls, indices);
    AddAmpls(ampls,indices);
    ampls->clear(); indices->clear();
  }
  delete ampls; ampls=NULL; delete indices; indices=NULL;
//   // get amplitude tensor
//   if(p_currents[0] && p_currents[1]) {                         // new for currents
//     ContractCurrents(p_currents[0], p_currents[1],
//                          p_momenta,Spin_Correlation_Tensor::Get_k0_n(),m_mefactor,
//                          p_ampls,p_indices);
//   }
//   else {
//     (*p_me)(  p_momenta,                                          // phase space point
//               p_ampls, p_indices,                                 // ampl. tensor and indices
//               Spin_Correlation_Tensor::Get_k0_n() );              // spinor base
//   }
  double value;
  if( p_ampls->size()==0 ) {                                    // no ampls <=> isotropic
    CreateTrivial(sigma);                                       // create trivial amplitude tensor
    value = 1.;
  }
  else {                                                        // not isotropic
    if( p_indices->size() ) {                                   // if there are indices
      Spin_Correlation_Tensor help_sct ( p_indices, p_ampls );  // temporary SCT
      help_sct.Contract(0,sigma);                               // contract over mother sigma
      value = real( help_sct.Trace() );                         // get T by taking trace of rest
    }
    else {                                                      // no spin correlation
      value = norm( (*p_ampls)[0] );                            // only one amplitude
    }
  }
  return value*weight*m_symmetry;
}

void Hadron_Decay_Channel::CreateTrivial( Spin_Density_Matrix * sigma )
{
  // create trivial amplitude tensor and its index bookkeeping
  p_indices->clear();
  p_ampls->clear();
  if( sigma ) {                                               // if spin > 0
    p_indices->push_back(pair<int,int>(0,sigma->Spin()));     // order of ind. does not matter
    for( int j=0; j<sigma->NrEntries(); ++j ) 
      p_ampls->push_back(1.);                                 // M=1 <=> isotropic
  }
  int spin;
  for( int i=0; i<m_nout; ++i ) {
    spin = p_flavours[i+1].IntSpin();                         // 2*spin of daughter
    if( spin ) {                                              // if spin > 0
      p_indices->push_back(pair<int,int>(i+1,spin));          // order of ind. does not matter
      for( int j=0; j<sqr(spin+1); ++j ) 
        p_ampls->push_back(1.);                               // M=1 <=> isotropic
    }
  }
}

void Hadron_Decay_Channel::AddAmpls( vector<Complex>* newampls, vector<pair<int,int> >* newindices )
{
    msg.Error()<<METHOD<<": Error. I am not able to handle multiple ME/current structures yet. "
      <<"FS is working on it though. Will abort..."<<endl;
    abort();
}


Current_Base* Hadron_Decay_Channel::SelectCurrent(string current_string)
{
// //   Heavy_Vector_Current[0,1]
  Data_Reader reader("[","]","#");
  reader.AddIgnore(",");
  vector<string> resultstrings;
  reader.VectorFromString(resultstrings,"",current_string,vtc::horizontal);
  
  Flavour_Info fi;
  fi.nout=resultstrings.size()-1;
  int* indices = new int[fi.nout];
  for(int i=0; i<fi.nout; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  fi.flavs=p_flavours; fi.indices=indices;

  Current_Base* current = Current_Getter_Function::GetObject(resultstrings[0],fi);
  if(current==NULL) {
  msg.Error()<<METHOD<<": Error. Current \""<<resultstrings[0]<<"\" specified in "
      <<m_path<<m_filename<<" was not recognized as a valid current. Will abort."<<endl;
    abort();
  }
  return current;
}


HD_ME_Base * Hadron_Decay_Channel::SelectME(string me_string)
{
  Flavour_Info fi;
  fi.nout=m_nout;
  fi.flavs=p_flavours;
  
  HD_ME_Base* me = HD_ME_Getter_Function::GetObject(me_string,fi);
  if(me==NULL) {
    msg.Error()<<METHOD<<": Error. Matrix element \""<<me_string<<"\" specified in "
      <<m_path<<m_filename<<" was not recognized as a valid ME. Will abort."<<endl;
    abort();
  }
  return me;
}

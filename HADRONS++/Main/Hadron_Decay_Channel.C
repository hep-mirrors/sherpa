#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "HADRONS++/ME_Library/Generic.H"
#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "ATOOLS/Phys/Decay_Table.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel( Flavour flin, string _path ) :
  Decay_Channel(flin),
  m_path(_path),
  p_amplitudes(NULL)
{
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
  if (p_ps) { delete p_ps; p_ps=NULL; }
  if(p_amplitudes) delete p_amplitudes; p_amplitudes=NULL;
  if(p_momenta)    delete [] p_momenta; p_momenta=NULL;
  if(p_flavours)   delete [] p_flavours; p_flavours=NULL;
  if(p_physicalflavours) delete [] p_physicalflavours; p_physicalflavours=NULL;
  for(size_t i=0;i<m_mes.size();i++)      delete m_mes[i].second;
  for(size_t i=0;i<m_currents.size();i++) {
    delete m_currents[i].second.first;
    delete m_currents[i].second.second;
  }
}


void Hadron_Decay_Channel::SetFileName(std::string filename)
{
  if(filename=="") {
    filename += GetDecaying().ShellName() + string("_");
    for ( int i=0; i<NOut(); i++ ) {
      filename += GetDecayProduct(i).ShellName();
    }
    filename += string(".dat");
  }
  m_filename = filename;
}


void Hadron_Decay_Channel::Initialise(GeneralModel startmd)
{
  p_flavours = new Flavour[GetN()];
  p_momenta  = new Vec4D[GetN()];
  p_flavours[0] = GetDecaying();
  int i=0;
  double totalmass=0.0;
  for (FlSetConstIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) {
    p_flavours[i+1] = (*flit);
    totalmass+=p_flavours[i+1].HadMass();
    i++;
  }
  if(totalmass>p_flavours[0].HadMass()) {
    msg_Error()<<"Error in "<<METHOD<<" for "<<Name()<<endl;
    msg_Error()<<"  Total outgoing mass heavier than incoming particle. Will abort."<<endl;
    abort();
  }
  p_physicalflavours=new Flavour[GetN()];
  for (int i=0; i<GetN(); ++i) {
    map<kf_code,kf_code>::const_iterator it =
        startmd.m_aliases.find(p_flavours[i].Kfcode());
    if (it!=startmd.m_aliases.end())
      p_physicalflavours[i] = Flavour(it->second, p_flavours[i].IsAnti());
    else
      p_physicalflavours[i] = p_flavours[i];
  }
  p_amplitudes = new Spin_Amplitudes(p_physicalflavours,GetN(),Complex(0.0,0.0));
  p_ps = new HD_PS_Base(this);
  // check for identical particles
  Flavour refflav;
  double symfactor(1.0);
  std::map<Flavour,size_t> fc;
  for (int i=1; i<1+NOut(); ++i) {
    std::map<Flavour,size_t>::iterator fit(fc.find(p_physicalflavours[i]));
    if (fit==fc.end()) fit=fc.insert(make_pair(p_physicalflavours[i],0)).first;
    ++fit->second;
  }
  for (std::map<Flavour,size_t>::const_iterator fit(fc.begin());
       fit!=fc.end();++fit) {
    symfactor*=Factorial(fit->second);
  }
  m_symmetry = 1./symfactor;
  m_cp_asymmetry_C=0.0; m_cp_asymmetry_S=0.0;

  m_startmd=startmd;

  // check if dc file exists
  My_In_File dcf(m_path,m_filename);
  if (dcf.Open()) {
    dcf.Close();
    msg_Tracking()<<METHOD<<": read "<<m_path<<m_filename<<endl;
    Data_Reader reader(" ",";","!");
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
      msg_Error()<<METHOD<<": Error.\n"
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
      msg_Error()<<METHOD<<": Error.\n"
                 <<"  Read in failure for <ME> section in "<<m_path
                 <<m_filename<<", will abort."<<endl;
      abort();
    }

    // process <Phasespace>
    vector<vector<string> > ps_svv;
    reader.SetFileBegin("<Phasespace>"); reader.SetFileEnd("</Phasespace>");
    reader.RereadInFile();
    if(!reader.MatrixFromFile(ps_svv)) {
    msg_Error()<<METHOD<<": Error."<<"  Read in failure for <Phasespace> section in "
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
    int* decayindices = new int[NOut()+1];
    for(int i=0;i<NOut()+1;i++) {
      decayindices[i]=i;
    }
    m_mes.push_back(MEPair(1.0,new Generic(p_physicalflavours,NOut()+1,
                                           decayindices,"Generic")));
    delete [] decayindices;
    p_ps->AddChannel( string("Isotropic"), 1., m_startmd );
    vector<double> results = CalculateResults();
    SetAlwaysIntegrate(false);
    WriteOut(results, true);
  }
}

void Hadron_Decay_Channel::ProcessOptions(vector<vector<string> > helpsvv)
{
  for (size_t i=0;i<helpsvv.size();i++) {
    if (helpsvv[i][0]==string("AlwaysIntegrate")) {
      SetAlwaysIntegrate( atoi( helpsvv[i][2].c_str() ) );
    }
    else if (helpsvv[i][0]==string("CPAsymmetryS")) {
      m_cp_asymmetry_S = ToType<double>(helpsvv[i][2]);
    }
    else if (helpsvv[i][0]==string("CPAsymmetryC")) {
      m_cp_asymmetry_C = ToType<double>(helpsvv[i][2]);
    }
  }
  // convert C and S to lambda, assuming DeltaGamma=0 for the determination of C and S.
  // this allows DeltaGamma dependent terms in the asymmetry
  double Abs2 = -1.0*(m_cp_asymmetry_C-1.0)/(m_cp_asymmetry_C+1.0);
  double Im = m_cp_asymmetry_S/(m_cp_asymmetry_C+1.0);
  double Re = sqrt(Abs2-sqr(Im));
  m_cp_asymmetry_lambda = Complex(Re, Im);
}

void Hadron_Decay_Channel::ProcessPhasespace(vector<vector<string> > ps_svv,
                                             Data_Reader           & reader,
                                             GeneralModel const    & model_for_ps)
{
  int nr_of_channels=0;
  for (size_t i=0;i<ps_svv.size();i++) {
    double weight = ToType<double>(ps_svv[i][0]);
    if( p_ps->AddChannel( ps_svv[i][1], weight, model_for_ps ) ) nr_of_channels++;
    else {
      msg_Error()<<METHOD<<": Warning. "<<ps_svv[i][1]<<" in "<<m_path<<m_filename
        <<" is not a valid phase space channel. Will ignore it."<<endl;
    }
  }
  if(nr_of_channels == 0) {
    msg_Error()<<METHOD<<": Warning. No valid phase space channels found in "
      <<m_path<<m_filename<<". Using Isotropic."<<endl;
    p_ps->AddChannel( string("Isotropic"), 1., m_startmd );
  }
}

void Hadron_Decay_Channel::ProcessME( vector<vector<string> > me_svv,
                                      Data_Reader           & reader,
                                      GeneralModel          & model_for_ps )
{
  int nr_of_mes=0;
  Algebra_Interpreter ip;
  ip.AddTag("GF", "8.24748e-6");
  for (size_t i=0;i<me_svv.size();i++) {
    if(me_svv[i].size()==3) {
      msg_Tracking()<<"Selecting ME for "<<Name()<<endl;
      HD_ME_Base* me = SelectME( me_svv[i][2] );
      me->SetPath(m_path);
      msg_Tracking()<<"  "<<me->Name()<<endl;
      vector<vector<string> > parameter_svv;
      reader.SetFileBegin("<"+me_svv[i][2]+">"); reader.SetFileEnd("</"+me_svv[i][2]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter_svv);
      GeneralModel me_model = Parameters2Model(parameter_svv, &model_for_ps);
      me->SetModelParameters( me_model );
      Complex factor = Complex(ToType<double>(ip.Interprete(me_svv[i][0])),
                               ToType<double>(ip.Interprete(me_svv[i][1])));
      m_mes.push_back( MEPair(factor,me) );
      nr_of_mes++;
    }
    if(me_svv[i].size()==4) {
      msg_Tracking()<<"Selecting currents for "<<Name()<<endl;
      Current_Base* current1 = SelectCurrent(me_svv[i][2]);
      current1->SetPath(m_path);
      vector<vector<string> > parameter1_svv;
      reader.SetFileBegin("<"+me_svv[i][2]+">"); reader.SetFileEnd("</"+me_svv[i][2]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter1_svv);
      GeneralModel current1_model = Parameters2Model(parameter1_svv, &model_for_ps);
      current1->SetModelParameters( current1_model );

      Current_Base* current2 = SelectCurrent(me_svv[i][3]);
      current2->SetPath(m_path);
      vector<vector<string> > parameter2_svv;
      reader.SetFileBegin("<"+me_svv[i][3]+">"); reader.SetFileEnd("</"+me_svv[i][3]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter2_svv);
      GeneralModel current2_model = Parameters2Model(parameter2_svv, &model_for_ps);
      current2->SetModelParameters( current2_model );

      msg_Tracking()<<"  "<<current1->Name()<<endl;
      msg_Tracking()<<"  "<<current2->Name()<<endl;

      // Sanity checks for current selection
      if( int(GetN()) != current1->GetN() + current2->GetN() ) {
        msg_Error()<<"Error in "<<METHOD<<": Current selection does not look sane "
                   <<"for "<<Name()<<". Check decaychannelfile."<<std::endl;
        abort();
      }

      Complex factor = Complex(ToType<double>(ip.Interprete(me_svv[i][0])),
                               ToType<double>(ip.Interprete(me_svv[i][1])));

      m_currents.push_back( CurrentsPair( factor, make_pair(current1,current2) ) );
      nr_of_mes++;
    }
  }
  if(nr_of_mes == 0) {
  msg_Error()<<METHOD<<": Warning. No valid matrix element found in "
      <<m_path<<m_filename<<". Using Generic."<<endl;
    int* decayindices = new int[NOut()+1];
    for(int i=0;i<NOut()+1;i++) {
      decayindices[i]=i;
    }
    m_mes.push_back(MEPair(1.0,new Generic(p_physicalflavours,NOut()+1,
                                           decayindices,"Generic")));
    delete [] decayindices;
  }
}

void Hadron_Decay_Channel::ProcessResult(vector<vector<string> > result_svv)
{
  if(result_svv.size()!=1) {
    vector<double> results = CalculateResults();
    WriteOut(results);
  }
  else if(result_svv[0].size()==3 && m_always_integrate) {
    vector<double> results = CalculateResults();
    // check whether result is different from before and write out if it is
    double oldwidth=ToType<double>(result_svv[0][0]);
    double oldmax=ToType<double>(result_svv[0][2]);
    if(oldwidth!=results[0] || oldmax!=results[2])
      WriteOut(results);
  }
  else if(result_svv[0].size()==3) {
    p_ps->SetResult(ToType<double>(result_svv[0][0]));
    p_ps->SetError(ToType<double>(result_svv[0][1]));
    p_ps->SetMaximum(ToType<double>(result_svv[0][2]));
  }
  else
    THROW(fatal_error, "Result section of "+m_path+"/"+m_filename+" did not "+
          "contain three entries. Aborting.");
}

GeneralModel Hadron_Decay_Channel::Parameters2Model(vector<vector<string> > helpsvv,
                                                    GeneralModel          * other_model )
{
  GeneralModel model(m_startmd);
  Algebra_Interpreter ip;
  for (size_t i=0;i<helpsvv.size();i++) {
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

bool Hadron_Decay_Channel::WriteOut( vector<double> results, bool newfile ) {

  if ( newfile ) {                // if DC file doesn't exist yet
    ofstream to;
    to.open((m_path+m_filename).c_str(),ios::out);

    // write header
    to<<"# Decay: "<<Name()<<endl;
    to<<"#        "<<setw(m_flin.IDName().length())<<left<<"0"<<" --> ";
    int i=0;
    for (FlSetConstIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) {
      to<<setw(flit->IDName().length()+1)<<left<<i+1;
      i++;
    }
    to<<endl<<endl;

    // write out options
    to<<"<Options>"<<endl;
    to<<"  AlwaysIntegrate = "<<m_always_integrate<<"    # 0...read results and skip integration"<<endl;
    to<<"                         # 1...don't read results and integrate"<<endl;
    to<<"</Options>"<<endl<<endl;

    // write out phasespace settings
    to<<"<Phasespace>"<<endl;
    to<<"  1.0 Isotropic"<<endl;
    to<<"</Phasespace>"<<endl<<endl;

    // write out ME settings
    to<<"<ME>"<<endl;
    to<<"  1.0 0.0 Generic"<<endl;
    to<<"</ME>"<<endl<<endl;

    // write out result
    to<<"<Result>"<<endl;
    int oldprec=to.precision(4);
    to<<"  "<<results[0]<<" "<<results[1]<<" "<<results[2]<<";"<<endl;
    to.precision(oldprec);
    to<<"</Result>"<<endl;
    to.close();
  } // if (read DC file)
  else {                                // if DC file exists
    Move(m_path+m_filename, m_path+"."+m_filename+".old");
    ofstream to((m_path+m_filename).c_str(),ios::out);

    // copy Options, Phasespace, ME, ...
    char buffer[100];
    ifstream from;
    from.open((m_path+"."+m_filename+string(".old")).c_str());
    bool extra_line = true;
    while (from.getline(buffer,100)) {
      if (buffer==string("<Result>")) { extra_line=false; break; }
      else to<<buffer<<endl;
    }
    from.close();

    // write out result
    if(extra_line) to<<endl;
    to<<"<Result>"<<endl;
    int oldprec=to.precision(4);
    to<<"  "<<results[0]<<" "<<results[1]<<" "<<results[2]<<";"<<endl;
    to.precision(oldprec);
    to<<"</Result>"<<endl;
    to.close();
  }
  return 1;
}

void Hadron_Decay_Channel::CalculateAmplitudes(Vec4D* moms, 
                                               Spin_Amplitudes* amps, bool anti)
{
  for(size_t i(0); i<amps->size(); ++i) (*amps)[i]=Complex(0.0,0.0);
  for(vector<MEPair>::iterator mpit=m_mes.begin(); mpit!=m_mes.end();mpit++) {
    mpit->second->SetAnti(anti);
    (*mpit->second)(moms, amps);
    for (size_t i(0);i<amps->size();++i)
      (*amps)[i]=(*amps)[i]*mpit->first;
  }

  for(vector<CurrentsPair>::iterator cpit=m_currents.begin(); 
      cpit!=m_currents.end();cpit++) {
    (*cpit).second.first->SetAnti(anti);
    (*cpit).second.second->SetAnti(anti);
    ContractCurrents((*cpit).second.first, (*cpit).second.second,
                     moms, (*cpit).first, amps);
  }
}

void Hadron_Decay_Channel::CalculateAmplitudes(Vec4D* moms, Amplitude_Tensor* amps, bool anti)
{
  Spin_Amplitudes* spinamps = new Spin_Amplitudes(amps->Particles());
  CalculateAmplitudes(moms, spinamps, anti);
  for(size_t i=0;i<amps->size();i++) amps->Insert(std::vector<Complex>(1,(*spinamps)[i]),i);
  delete spinamps;
}

vector<double> Hadron_Decay_Channel::CalculateResults()
{
  long int seed = ran->GetSeed();
  ran->SetSeed(123456);
  vector<double> results = p_ps->CalculateNormalisedWidth();
  ran->SetSeed(seed);
  return results;
}

bool Hadron_Decay_Channel::SetColorFlow(Particle_Vector outparts,int n_q, int n_g)
{
  if(n_q==0 && n_g==0)  return true;
  if(m_mes.size()>0) {
    // try if the matrix element knows how to set the color flow
    if(m_mes[0].second->SetColorFlow(outparts,n_q,n_g)) return true;
  }
  // otherwise try some common situations
  int n=outparts.size();
  if(n_q==2 && n_g==0 && n==2) {
    if(outparts[0]->Flav().IsAnti()) {
      outparts[0]->SetFlow(2,-1);
      outparts[1]->SetFlow(1,outparts[0]->GetFlow(2));
    }
    else {
      outparts[0]->SetFlow(1,-1);
      outparts[1]->SetFlow(2,outparts[0]->GetFlow(1));
    }
    return true;
  }
  else if(n_q==0 && n_g==2) {
    int inflow(-1), outflow(-1);
    Particle_Vector::iterator pit;
    for(pit=outparts.begin(); pit!=outparts.end(); pit++) {
      if((*pit)->Flav().IsGluon()) {
        (*pit)->SetFlow(2,inflow);
        (*pit)->SetFlow(1,outflow);
        inflow=(*pit)->GetFlow(1);
        outflow=(*pit)->GetFlow(2);
      }
    }
    return true;
  }
  else if(n_q==0 && n_g==n) {
    outparts[0]->SetFlow(2,-1);
    outparts[0]->SetFlow(1,-1);
    for(int i=1;i<n-1;++i) {
      unsigned int c=Flow::Counter();
      outparts[i]->SetFlow(2,c-1);
      outparts[i]->SetFlow(1,c);
    }
    outparts[n-1]->SetFlow(2,outparts[n-2]->GetFlow(1));
    outparts[n-1]->SetFlow(1,outparts[0]->GetFlow(2));
    return true;
  }
  else return false;
}

// differential with random PS points; just for weight
double Hadron_Decay_Channel::Differential()
{
  p_momenta[0] = Vec4D(p_physicalflavours[0].HadMass(),0.,0.,0.);        // decay from rest
  p_ps->GeneratePoint(p_momenta,(PHASIC::Cut_Data*)(NULL));     // generate a PS point
  p_ps->GenerateWeight(p_momenta,(PHASIC::Cut_Data*)(NULL));    // calculate its weight factor
  double weight = p_ps->Weight();                               // get weight factor
  CalculateAmplitudes(p_momenta,p_amplitudes,false);
  double value=p_amplitudes->SumSquare();
  value /= (GetDecaying().IntSpin()+1);
  value *= m_symmetry;
  return value*weight;
}

// differential with incoming momentum; for weight and momenta
double Hadron_Decay_Channel::Differential(Vec4D * mom, bool anti)
{
#ifdef DEBUG__Hadrons
  if( !IsZero(mom[0][1]) || !IsZero(mom[0][2]) || !IsZero(mom[0][3]) ) {
    PRINT_INFO("Error: given momentum is not in CMS: "<<mom[0]);
  }
#endif
  p_ps->GeneratePoint(mom,(PHASIC::Cut_Data*)(NULL));     // generate a PS point
  p_ps->GenerateWeight(mom,(PHASIC::Cut_Data*)(NULL));    // calculate its weight factor
  
  double weight = p_ps->Weight();                               // weight factor
  CalculateAmplitudes(mom,p_amplitudes,anti);
  double value=p_amplitudes->SumSquare();
  value /= (GetDecaying().IntSpin()+1);
  return value*weight;
}


Current_Base* Hadron_Decay_Channel::SelectCurrent(string current_string)
{
  Data_Reader reader(",",";","#","]");
  reader.AddWordSeparator("[");
  vector<string> resultstrings;
  reader.SetString(current_string);
  reader.VectorFromString(resultstrings);
  
  Flavour_Info fi;
  fi.n=resultstrings.size()-1;
  int* indices = new int[fi.n]; // will get deleted by Current_Base
  for(int i=0; i<fi.n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  fi.flavs=p_physicalflavours; fi.indices=indices;

  Current_Base* current = Current_Getter_Function::GetObject(resultstrings[0],fi);
  if(current==NULL) {
  msg_Error()<<METHOD<<": Error. Current \""<<resultstrings[0]<<"\" specified in "
      <<m_path<<m_filename<<" was not recognized as a valid current. Will abort."<<endl;
    abort();
  }
  return current;
}


HD_ME_Base * Hadron_Decay_Channel::SelectME(string me_string)
{
  Data_Reader reader(",",";","#","]");
  reader.AddWordSeparator("[");
  vector<string> resultstrings;
  reader.SetString(me_string);
  reader.VectorFromString(resultstrings);
  if(resultstrings.size()==1 && resultstrings[0]=="Generic") {
    for(int i=0;i<NOut()+1;i++) 
      resultstrings.push_back( ToString<size_t>(i) );
  }
  if(int(resultstrings.size())!=NOut()+2) {
    msg_Error()<<METHOD<<" Error: Number of indices in \""<<me_string<<"\" ("
      <<int(resultstrings.size())-1<<") in "<<m_path<<m_filename<<" doesn't "
      <<"equal number of particles ("<<NOut()+1<<"). Will abort."<<endl;
    abort();
  }

  Flavour_Info fi;
  fi.n=NOut()+1;
  int* indices = new int[fi.n]; // will get deleted by HD_ME_Base
  for(int i=0; i<fi.n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  fi.flavs=p_physicalflavours; fi.indices=indices;

  HD_ME_Base* me = HD_ME_Getter_Function::GetObject(resultstrings[0],fi);
  if(me==NULL) {
    msg_Error()<<METHOD<<": Error. Matrix element \""<<me_string<<"\" specified in "
      <<m_path<<m_filename<<" was not recognized as a valid ME. Will abort."<<endl;
    abort();
  }
  return me;
}


namespace ATOOLS {
  template <> Blob_Data<HADRONS::Hadron_Decay_Channel*>::~Blob_Data() { }
  template class Blob_Data<HADRONS::Hadron_Decay_Channel*>;
  template HADRONS::Hadron_Decay_Channel* &Blob_Data_Base::Get<HADRONS::Hadron_Decay_Channel*>();
}


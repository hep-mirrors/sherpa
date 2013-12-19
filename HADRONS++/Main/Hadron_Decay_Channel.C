#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "HADRONS++/ME_Library/Generic.H"
#include "HADRONS++/ME_Library/Current_ME.H"
#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Color_Function_Decay.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Blob.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel(Flavour fl, const Mass_Selector* ms,
                                           string _path) :
  Decay_Channel(fl, ms),
  m_path(_path), m_always_integrate(false),
  m_cp_asymmetry_C(0.0), m_cp_asymmetry_S(0.0)
{
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
}


void Hadron_Decay_Channel::SetFileName(std::string filename)
{
  if(filename=="") {
    filename += GetDecaying().ShellName();
    if (filename=="B_{s}") filename = "Bs";
    filename += string("_");
    for ( int i=0; i<NOut(); i++ ) {
      filename += GetDecayProduct(i).ShellName();
    }
    filename += string(".dat");
  }
  m_filename = filename;
}


void Hadron_Decay_Channel::Initialise(GeneralModel startmd)
{
  m_physicalflavours=m_flavours;
  for (size_t i=0; i<m_flavours.size(); ++i) {
    map<kf_code,kf_code>::const_iterator it = 
        startmd.m_aliases.find(m_flavours[i].Kfcode());
    if (it!=startmd.m_aliases.end())
      m_physicalflavours[i] = Flavour(it->second, m_flavours[i].IsAnti());
  }
  
  double totalmass=0.0;
  for (size_t i=1; i<m_flavours.size(); ++i) {
    totalmass+=m_flavours[i].HadMass();
  }
  if(totalmass>m_flavours[0].HadMass()) {
    msg_Error()<<"Error in "<<METHOD<<" for "<<Name()<<"\n"
	       <<"    Total outgoing mass heavier than incoming particle.\n"
	       <<"    Will abort.\n";
    abort();
  }
  SetChannels(new PHASIC::Multi_Channel(""));
  Channels()->SetNin(1);
  Channels()->SetNout(NOut());
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
		 <<"   Read in failure for <Options> section in "
		 <<m_path<<m_filename<<".\n"
		 <<"   Will abort."<<endl;
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
    msg_Error()<<METHOD<<": Error.\n"
	       <<"   Read in failure for <Phasespace> section in "
	       <<m_path<<m_filename<<".\n"
	       <<"   Will abort."<<endl;
      abort();
    }
    ProcessPhasespace(ps_svv, reader, model_for_ps);

    // process <Result> 
    // don't do it before ME and phasespace, or CalcNormWidth doesn't work!
    vector<vector<string> > result_svv;
    reader.SetFileBegin("<Result>"); reader.SetFileEnd("</Result>");
    reader.RereadInFile();
    reader.MatrixFromFile(result_svv);
    ProcessResult(result_svv);
  }
  else { // if DC file does not exist yet
    msg_Tracking()<<"No DC file yet in :"<<m_path<<"/"<<m_filename<<".\n";
    int n=NOut()+1;
    vector<int> decayindices(n);
    for(int i=0;i<n;i++) decayindices[i]=i;
    HD_ME_Base* me=new Generic(m_physicalflavours,decayindices,"Generic");
    PHASIC::Color_Function_Decay* col=new PHASIC::Color_Function_Decay();
    AddDiagram(me, col);
    AddPSChannel( string("Isotropic"), 1., m_startmd);
    msg_Tracking()<<"Calculating width for "<<Name()<<":\n";
    CalculateWidth();
    msg_Tracking()<<"   yields "<<m_iwidth<<".\n";
    WriteOut(true);
  }
}

void Hadron_Decay_Channel::ProcessOptions(vector<vector<string> > helpsvv)
{
  for (size_t i=0;i<helpsvv.size();i++) {
    if (helpsvv[i][0]==string("AlwaysIntegrate")) {
      m_always_integrate=atoi(helpsvv[i][2].c_str());
    }
    else if (helpsvv[i][0]==string("CPAsymmetryS")) {
      m_cp_asymmetry_S = ToType<double>(helpsvv[i][2]);
    }
    else if (helpsvv[i][0]==string("CPAsymmetryC")) {
      m_cp_asymmetry_C = ToType<double>(helpsvv[i][2]);
    }
  }
  // convert C and S to lambda, assuming DeltaGamma=0 for the determination 
  // of C and S.
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
    if(AddPSChannel( ps_svv[i][1], weight, model_for_ps ) ) nr_of_channels++;
    else {
      msg_Error()<<METHOD<<":  Warning\n"
		 <<"   "<<ps_svv[i][1]<<" in "<<m_path<<m_filename
		 <<" is not a valid phase space channel.\n"
		 <<"   Will ignore it and hope for the best.\n";
    }
  }
  if(nr_of_channels == 0) {
    msg_Error()<<METHOD<<": Warning. No valid phase space channels found in "
	       <<m_path<<m_filename<<". Using Isotropic."<<endl;
    AddPSChannel( string("Isotropic"), 1., m_startmd );
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
      GeneralModel me_model=Parameters2Model(parameter_svv,model_for_ps);
      me->SetModelParameters( me_model );
      Complex factor = Complex(ToType<double>(ip.Interprete(me_svv[i][0])),
                               ToType<double>(ip.Interprete(me_svv[i][1])));
      me->SetFactor(factor);
      PHASIC::Color_Function_Decay* col=new PHASIC::Color_Function_Decay();
      AddDiagram(me, col);
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
      GeneralModel current1_model=Parameters2Model(parameter1_svv,model_for_ps);
      current1->SetModelParameters( current1_model );

      Current_Base* current2 = SelectCurrent(me_svv[i][3]);
      current2->SetPath(m_path);
      vector<vector<string> > parameter2_svv;
      reader.SetFileBegin("<"+me_svv[i][3]+">"); reader.SetFileEnd("</"+me_svv[i][3]+">");
      reader.RereadInFile();
      reader.MatrixFromFile(parameter2_svv);
      GeneralModel current2_model=Parameters2Model(parameter2_svv,model_for_ps);
      current2->SetModelParameters( current2_model );

      msg_Tracking()<<"  "<<current1->Name()<<endl;
      msg_Tracking()<<"  "<<current2->Name()<<endl;

      // Sanity checks for current selection
      if(size_t(1+NOut()) != current1->DecayIndices().size()+
         current2->DecayIndices().size()) {
        msg_Error()<<"Error in "<<METHOD<<": Current selection does not look sane "
                   <<"for "<<Name()<<". Check decaychannelfile."<<std::endl;
        abort();
      }

      Complex factor = Complex(ToType<double>(ip.Interprete(me_svv[i][0])),
                               ToType<double>(ip.Interprete(me_svv[i][1])));

      vector<int> indices (NOut()+1);
      for(int i=0; i<NOut()+1; i++) indices[i] = i;

      Current_ME* me=
        new Current_ME(m_physicalflavours, indices, "Current_ME");
      me->SetCurrent1(current1);
      me->SetCurrent2(current2);
      me->SetFactor(factor);
      PHASIC::Color_Function_Decay* col=new PHASIC::Color_Function_Decay();
      AddDiagram(me, col);
      nr_of_mes++;
    }
  }
  if(nr_of_mes == 0) {
    msg_Error()<<METHOD<<": Warning. No valid matrix element found in "
               <<m_path<<m_filename<<". Using Generic."<<endl;
    int n=NOut()+1;
    vector<int> decayindices(n);
    for(int i=0;i<n;i++) decayindices[i]=i;
    HD_ME_Base* me=new Generic(m_physicalflavours,decayindices,"Generic");
    PHASIC::Color_Function_Decay* col=new PHASIC::Color_Function_Decay();
    AddDiagram(me, col);
  }
}

void Hadron_Decay_Channel::ProcessResult(vector<vector<string> > result_svv)
{
  if(result_svv.size()!=1) {
    msg_Info()<<"Calculating width (PR1) for "<<Name()<<endl;
    CalculateWidth();
    msg_Info()<<"   yields "<<m_iwidth<<".\n";
    WriteOut();
  }
  else if(result_svv[0].size()==3 && m_always_integrate) {
    msg_Info()<<"Calculating width (PR2) for "<<Name()<<endl;
    CalculateWidth();
    msg_Info()<<"   yields "<<m_iwidth<<".\n";
    // check whether result is different from before and write out if it is
    double oldwidth=ToType<double>(result_svv[0][0]);
    double oldmax=ToType<double>(result_svv[0][2]);
    if(oldwidth!=m_iwidth || oldmax!=m_max)
      WriteOut();
  }
  else if(result_svv[0].size()==3) {
    m_iwidth=ToType<double>(result_svv[0][0]);
    m_ideltawidth=ToType<double>(result_svv[0][1]);
    m_max=ToType<double>(result_svv[0][2]);
  }
  else
    THROW(fatal_error, "Result section of "+m_path+"/"+m_filename+" did not "+
          "contain three entries. Aborting.");
}

GeneralModel
Hadron_Decay_Channel::Parameters2Model(vector<vector<string> > helpsvv,
                                       GeneralModel& other_model)
{
  GeneralModel model(m_startmd);
  Algebra_Interpreter ip;
  for (size_t i=0;i<helpsvv.size();i++) {
    if ( helpsvv[i][1] == string("=")) {
      if( helpsvv[i].size() == 3 ) {        // <name> = <real value>
        double real = ToType<double> (ip.Interprete(helpsvv[i][2]) );
        model[helpsvv[i][0]] = real;
        other_model[helpsvv[i][0]] = real;
      }
      if( helpsvv[i].size() == 4 ) {        // <name> = <complex value>
        double abs   = ToType<double>(ip.Interprete(helpsvv[i][2]) );
        double phase = ToType<double>( ip.Interprete(helpsvv[i][3]) );
        model[helpsvv[i][0]+string("_abs")] = abs;
        model[helpsvv[i][0]+string("_phase")] = phase;
        other_model[helpsvv[i][0]+string("_abs")] = abs;
        other_model[helpsvv[i][0]+string("_phase")] = phase;
      }
    }
    if ( helpsvv[i][2] == string("=")) {
      if( helpsvv[i].size() == 4 ) {        // <name> <index> = <real value>
        double real = ToType<double>( ip.Interprete(helpsvv[i][3]) );
        model[helpsvv[i][0]+string("_")+helpsvv[i][1]] = real;
        other_model[helpsvv[i][0]+string("_")+helpsvv[i][1]] = real;
      }
      if( helpsvv[i].size() == 5 ) {        // <name> <index> = <complex value>
        double abs   = ToType<double> ( ip.Interprete(helpsvv[i][3]) );
        double phase = ToType<double> ( ip.Interprete(helpsvv[i][4]) );
        model[helpsvv[i][0]+"_"+helpsvv[i][1]+"_abs"] = abs;
        model[helpsvv[i][0]+"_"+helpsvv[i][1]+"_phase"] = phase;
        other_model[helpsvv[i][0]+"_"+helpsvv[i][1]+"_abs"]=abs;
        other_model[helpsvv[i][0]+"_"+helpsvv[i][1]+"_phase"]=phase;
      }
    }
  }
  return model;
}

void Hadron_Decay_Channel::WriteOut(bool newfile) {
  if ( newfile ) {                // if DC file doesn't exist yet
    My_Out_File to(m_path+m_filename);
    to.Open();

    // write header
    *to<<"# Decay: "<<Name()<<endl;
    *to<<"#        "<<setw(m_flavours[0].IDName().length())<<left<<"0"<<" --> ";
    int i=0;
    for (size_t i=1; i<m_flavours.size(); ++i) {
      *to<<setw(m_flavours[i].IDName().length()+1)<<left<<i+1;
      i++;
    }
    *to<<endl<<endl;

    // write out options
    *to<<"<Options>"<<endl;
    *to<<"  AlwaysIntegrate = "<<m_always_integrate
       <<"    # 0...read results and skip integration"<<endl;
    *to<<"                         # 1...don't read results and integrate"<<endl;
    *to<<"</Options>"<<endl<<endl;

    // write out phasespace settings
    *to<<"<Phasespace>"<<endl;
    *to<<"  1.0 Isotropic"<<endl;
    *to<<"</Phasespace>"<<endl<<endl;

    // write out ME settings
    *to<<"<ME>"<<endl;
    *to<<"  1.0 0.0 Generic"<<endl;
    *to<<"</ME>"<<endl<<endl;

    // write out result
    *to<<"<Result>"<<endl;
    int oldprec=to->precision(4);
    *to<<"  "<<m_iwidth<<" "<<m_ideltawidth<<" "<<m_max<<";"<<endl;
    to->precision(oldprec);
    *to<<"</Result>"<<endl;
    to.Close();
  } // if (read DC file)
  else {                                
    // if DC file exists
    PRINT_INFO("TODO: migrate to Decaydata.db for "<<m_path+m_filename<<": "<<m_iwidth<<" "<<m_ideltawidth<<" "<<m_max);
    /* TODO migrate to Decaydata.db
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
    to<<"  "<<m_iwidth<<" "<<m_ideltawidth<<" "<<m_max<<";"<<endl;
    to.precision(oldprec);
    to<<"</Result>"<<endl;
    to.close();
    */
  }
}

bool Hadron_Decay_Channel::SetColorFlow(ATOOLS::Blob* blob)
{
  int n_q(0), n_g(0);
  for(int i=0;i<blob->NOutP();i++) {
    if(blob->OutParticle(i)->Flav().IsQuark())      n_q++;
    else if(blob->OutParticle(i)->Flav().IsGluon()) n_g++;
  }
  if(n_q>0 || n_g>0) {
    blob->SetStatus(blob_status::needs_showers);
    Particle_Vector outparts=blob->GetOutParticles();
    if(m_diagrams.size()>0) {
      // try if the matrix element knows how to set the color flow
      HD_ME_Base* firstme=(HD_ME_Base*) m_diagrams[0].first;
      bool anti=blob->InParticle(0)->Flav().IsAnti();
      if(firstme->SetColorFlow(outparts,n_q,n_g,anti)) return true;
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
    else {
      msg_Error()<<METHOD<<" wasn't able to set the color flow for"<<endl
                 <<*blob<<endl;
      return false;
    }
  }
  else return true;
}


Current_Base* Hadron_Decay_Channel::SelectCurrent(string current_string)
{
  Data_Reader reader(",",";","#","]");
  reader.AddWordSeparator("[");
  vector<string> resultstrings;
  reader.SetString(current_string);
  reader.VectorFromString(resultstrings);
  
  int n=resultstrings.size()-1;
  vector<int> indices(n);
  for(int i=0; i<n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  ME_Parameters fi(m_physicalflavours, indices);

  Current_Base* current=Current_Getter_Function::GetObject(resultstrings[0],fi);
  if(current==NULL) {
    msg_Error()<<METHOD<<": Current '"<<resultstrings[0]<<"' specified in "
               <<m_path<<m_filename<<" was not recognized as a valid current. "
               <<"Will abort."<<endl;
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

  int n=NOut()+1;
  vector<int> indices(n);
  for(int i=0; i<n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  ME_Parameters fi(m_physicalflavours, indices);

  HD_ME_Base* me = HD_ME_Getter_Function::GetObject(resultstrings[0],fi);
  if(me==NULL) {
    msg_Error()<<METHOD<<": Error. Matrix element \""<<me_string<<"\" specified in "
      <<m_path<<m_filename<<" was not recognized as a valid ME. Will abort."<<endl;
    abort();
  }
  return me;
}

void Hadron_Decay_Channel::LatexOutput(std::ostream& f, double totalwidth)
{
  f<<"$"<<GetDecaying().TexName()<<"$ $\\to$ ";
  for (size_t i=1; i<m_flavours.size(); ++i)
    f<<"$"<<m_flavours[i].TexName()<<"$ ";
  f<<" & ";
  char helpstr[100];
  sprintf( helpstr, "%.4f", Width()/totalwidth*100. );
  f<<helpstr;
  if( DeltaWidth() > 0. ) {
    sprintf( helpstr, "%.4f", DeltaWidth()/totalwidth*100. );
    f<<" $\\pm$ "<<helpstr;
  }
  f<<" \\% ";
  if(Origin()!="") {
    f<<"[\\verb;"<<Origin()<<";]";
  }
  f<<"\\\\"<<endl;
  if((m_diagrams.size()>0 &&
      ((HD_ME_Base*) m_diagrams[0].first)->Name()!="Generic")) {
    sprintf( helpstr, "%.4f", IWidth()/totalwidth*100. );
    f<<" & "<<helpstr;
    if( IDeltaWidth() > 0. ) {
      sprintf( helpstr, "%.4f", IDeltaWidth()/totalwidth*100. );
      f<<" $\\pm$ "<<helpstr;
    }
    f<<" \\% ";
  }
  for(size_t i=0;i<m_diagrams.size();i++) {
    HD_ME_Base* me=(HD_ME_Base*) m_diagrams[i].first;
    if(me->Name()=="Current_ME") {
      Current_ME* cme=(Current_ME*) me;
      f<<"\\verb;"<<cme->GetCurrent1()->Name()
       <<";$\\otimes$\\verb;"<<cme->GetCurrent2()->Name()<<"; & \\\\"<<endl;
    }
    else if (me->Name()=="Generic") {
      // do nothing
    }
    else {
      f<<"\\verb;"<<me->Name()<<"; & \\\\"<<endl;
    }
  }
}

bool Hadron_Decay_Channel::AddPSChannel(string name,double weight,
                                        GeneralModel const & md)
{
  PHASIC::Single_Channel * sc=
    HD_Channel_Selector::GetChannel(1, NOut(),&m_flavours.front(),name,md,p_ms);
  if (sc!=NULL) {
    sc->SetAlpha(weight);
    Channels()->Add(sc);
    return true;
  }
  else return false;
}


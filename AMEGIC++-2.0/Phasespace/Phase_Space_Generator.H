#ifndef Phase_Space_Generator_h
#define Phase_Space_Generator_h

#include "Multi_Channel.H"
#include "MyString.H"
#include "Point.H"

namespace AMEGIC {
  class Process_Base;

  class Phase_Space_Generator {
    static const int buffersize = 512*32;
    int         nin,nout;
    char        filename[12], procname[12];
    std::string path,pathID,pID;
    std::string headername,cfilename;
    int m_mode;

    bool  CheckForOldChannels(int,ATOOLS::Flavour *,PHASIC::Multi_Channel *);
    bool  IsFile(std::string &);
    bool  Search(std::ifstream &,std::string);
    int   Search(std::string file,std::string search);
    void  Copy(std::string,std::string);
    bool  MakeHeader(std::string &,char*,int);
    bool  MakeCfile(std::string &,char*,int,int);
    int   AddToHeader(std::string &,char*,bool);
    int   AddToCfile(std::string &,char*,int,bool);
    void  AddToMakefile(std::string Makefile,std::string pathID,std::string fileID);
    void  AddToMakefileAM(std::string, std::string, std::string);
    void  AddToSetChannel();
  public:
    Phase_Space_Generator(int,int);
    bool Construct(PHASIC::Multi_Channel *,std::string,std::string,ATOOLS::Flavour *, 
		   Process_Base *);
    PHASIC::Single_Channel * SetChannel(int,int,ATOOLS::Flavour* fl,int,std::string&);
  };

  //! The numbers of incoming and outgoing momenta.
  /*!
    Two identifiers : The path (a name for the subdirectory, like P{nin}_{nout}) and
    the specific process name, pID.
  */
  /*!
    Two more names for the generation of libraries. Somehow I have the impression, that
    we do some overcounting here, have to think about it.
  */
  /*!
    Obvious names for specific files.
  */
  /*!
    Checks wether the channels have already been constructed.
    If there's a mismatch of C and H files the program aborts.
    If the channel was already constructed, the corresponding Single_Channel is added 
    to the Multi_Channel and a 1 is returned, otherwise a zero is returned.
  */
  /*!
    Checks, if files with a given name are already existing. If this is the case, 1 is returned.
  */
  /*!
    Checks, if the string is contained in the stream. If so, a 1 is returned.
  */
  /*!
    Copies the ingredients of one file to another.
    The filenames are given by the strings.
  */
  /*!
    Create the H- and the C-file for the Multi_Channel.
    The basic idea for the Multi_Channel is that each Single_Channel
    receives a number according to which they will be called. 
    These numbers are stuffed into a switch-case structure
    that gets initialized in these methods.
  */
  /*!
    These methods add further channels into the structure of
    H-, C-, and Makefile as well as to the template SetChannel file. 
    The trick in all cases is to copy the old file, insert the
    lines to be included and copy back.
  */
  /*!
    This method constructs new channels for the phase space integrator.
    It checks, via ChechForOldChannels whether old channels are existing. 
    If this is the case Construct return a 0. Otherwise it proceeds with the 
    construction of channels with help of the Channel_Generator. If a channel was
    generated, i.e. some new files and methods were created, this has to be inserted 
    into the corresponding Makefiles, a link has to be established and the main headers 
    and cfiles (which distribute the individual channels) have to be modified. 
    At the end, if a new channel was created, "true" is returned.
  */
  /*!
    This method is located in an external file to be changed during the usage of AMEGIC. 
    It basically relates the library files for integration channels and the 
    Phase_Space_Handler/Phase_Space_Generator complex. Depending on the (name-) string
    and other chararcteristics it returns specific integration channels out of the
    channel library.
  */
}
#endif



#ifndef Phase_Space_Generator_h
#define Phase_Space_Generator_h

#include "Process_Base.H"
#include "Multi_Channel.H"
#include "ISR_Handler.H"
#include "MyString.H"
#include "Point.H"

namespace PHASIC {
  class Phase_Space_Generator {
    //! The numbers of incoming and outgoing momenta.
    int    nin,nout;
    /*!
      Two identifiers : The path (a name for the subdirectory, like P{nin}_{nout}) and
      the specific process name, pID.
    */
    std::string path,pathID,pID;
    /*!
      Two more names for the generation of libraries. Somehow I have the impression, that
      we do some overcounting here, have to think about it.
    */
    char   filename[12], procname[12];
    /*!
      Obvious names for specific files.
    */
    std::string headername,cfilename;
    static const int buffersize = 512*32;
    /*!
      Checks wether the channels have already been constructed.
      If there's a mismatch of C and H files the program aborts.
      If the channel was already constructed, the corresponding Single_Channel is added 
      to the Multi_Channel and a 1 is returned, otherwise a zero is returned.
    */
    bool  CheckForOldChannels(int,APHYTOOLS::Flavour *,Multi_Channel *);
    /*!
      Checks, if files with a given name are already existing. If this is the case, 1 is returned.
    */
    bool  IsFile(std::string &);
    /*!
      Checks, if the string is contained in the stream. If so, a 1 is returned.
    */
    bool  Search(std::ifstream &,std::string);
    int   Search(std::string file,std::string search);
    /*!
      Copies the ingredients of one file to another.
      The filenames are given by the strings.
    */
    void Copy(std::string,std::string);
    /*!
      Create the H- and the C-file for the Multi_Channel.
      The basic idea for the Multi_Channel is that each Single_Channel
      receives a number according to which they will be called. 
      These numbers are stuffed into a switch-case structure
      that gets initialized in these methods.
    */
    bool MakeHeader(std::string &,char*,int);
    bool MakeCfile(std::string &,char*,int,int);
    /*!
      These methods add further channels into the structure of
      H-, C-, and Makefile as well as to the template SetChannel file. 
      The trick in all cases is to copy the old file, insert the
      lines to be included and copy back.
    */
    int  AddToHeader(std::string &,char*,bool);
    int  AddToCfile(std::string &,char*,int,bool);
    void AddToMakefile(std::string Makefile,std::string pathID,std::string fileID);
  //    void AddToMakefile(std::string,std::string);
    void AddToSetChannel();
  public:
    Phase_Space_Generator(int,int);
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
    bool Construct(Multi_Channel *,std::string,std::string,APHYTOOLS::Flavour *, 
		   AMEGIC::Process_Base *);
    /*!
      This method creates the beam and the ISR channels for the phase space integration.
      The idea is to iterate over all type of the integer-vector and add
      in the corresponding channels. For Breit-Wigner type resonances, two more
      vectors are needed, carrying the mass and the width.
    */
    bool CreateBeamChannels(Multi_Channel *,APHYTOOLS::Flavour *,const std::vector<Channel_Info> &);
    bool CreateISRChannels(Multi_Channel *,APHYTOOLS::Flavour *,const std::vector<Channel_Info> &);
    bool CreateFSRChannels(Multi_Channel *,APHYTOOLS::Flavour *);
    /*!
      This method is located in an external file to be changed during the usage of AMEGIC. 
      It basically relates the library files for integration channels and the 
      Phase_Space_Handler/Phase_Space_Generator complex. Depending on the (name-) string
      and other chararcteristics it returns specific integration channels out of the
      channel library.
    */
    Single_Channel * SetChannel(int,int,APHYTOOLS::Flavour* fl,int,std::string&);
  };
}
#endif



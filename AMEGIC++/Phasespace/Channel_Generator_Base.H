#ifndef AMEGIC_Phasespace_Channel_Generator_Base_H
#define AMEGIC_Phasespace_Channel_Generator_Base_H

#include "PHASIC++/Channels/Single_Channel.H"
#include <string>
#include <vector>
#include <map>


namespace AMEGIC {

typedef std::map<std::string,std::string> Decls;
  class Point;
  
  typedef std::vector<std::string> String_List;

  class Channel_Generator_Base {  //: public PHASIC::Single_Channel {
  protected:
    int         nin,nout;
    Point     * plist;
    std::string GetMassIndex(std::string &str);
    std::string GetMassIndex(char &str);
  public:
    Channel_Generator_Base(int,int,Point *);
    virtual ~Channel_Generator_Base();
    virtual int  MakeChannel(int&,int,std::string&,std::string&) = 0;
    Point     * GetPointlist()             { return plist; }
    //void        SetName(std::string _name) { name = _name; }
    virtual std::string CreateChannelID(int) = 0;
    virtual int NumberOfChannels()         { return 1; }
  };
}  
#endif


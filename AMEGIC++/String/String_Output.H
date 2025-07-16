#ifndef AMEGIC_String_String_Output_H
#define AMEGIC_String_String_Output_H


#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/String/String_Tree.H"
#include "AMEGIC++/Main/Helicity.H"
#include "AMEGIC++/String/String_Library.H"

namespace AMEGIC {
class String_Output {
  String_Library slib;

  std::string path;
  int maxgraph,maxhel;
  std::string pID,pathID;

  int m_mode;

  void Cform(std::ofstream&,int,int,Virtual_String_Generator*,
	     sknot***,String_Tree*,const std::string&,Helicity*);
  void Zform(std::ofstream&,int,int,Virtual_String_Generator*,String_Tree*);
  void Make_Header(std::ofstream &,Virtual_String_Generator*);
  int  Line_Form(std::ofstream&,const std::string &);
  void Add_To_Set_Values();
public:
  String_Output(const std::string&,int,int,int); 
  void Output(sknot***,String_Tree*,Virtual_String_Generator*,Helicity*);
};
}
#endif




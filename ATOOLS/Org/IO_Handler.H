#ifndef ATOOLS_Org_IO_Handler_H
#define ATOOLS_Org_IO_Handler_H

#include <string>
#include <fstream>
#include <map>
#include <vector>

namespace ATOOLS {  

  typedef std::map<std::string,std::string> Variable_Map;
  typedef Variable_Map::iterator           Variable_Iterator;

  //! responsible for writing out whole classes
  class IO_Handler {
    std::fstream   m_file;
    std::string    m_filename;

    std::string    m_buffer;
    Variable_Map   m_vars;

    std::vector<char> m_seps, m_coms;

    unsigned long m_nx, m_ny;
  public:
    IO_Handler();

    ~IO_Handler();
    
    // set output filename
    int SetFileName(std::string);
    int SetFileNameRO(std::string);


    void SetSeparator(char);
    void AddSeparator(char);
    void SetComment(char);
    void AddComment(char);

    template <class Type> 
    int ValueInput(std::string name, Type &);
    
    // output file (compare rpa, etc.)
    template <class Type> 
    IO_Handler & operator<<(const Type &);

    // readin class from file 
    template <class Type> 
    IO_Handler & operator>>(Type &);

    //
    // standard output methods
    //

    template <class Type> 
    void MatrixOutput(const std::string name,Type ** const ,const int nx,const int ny);

    template <class Type> 
    void ArrayOutput(const std::string name,const Type *,const int nx,bool writesize=1);

    template <class Type> 
    void Output(const std::string name,const Type &);

    //
    // standard input methods
    //

    template <class Type> 
    Type ** MatrixInput(const std::string name=std::string(""), int nx=-1, int ny=-1);

    template <class Type> 
    Type * ArrayInput(const std::string name=std::string(""), int nx=-1);

    template <class Type> 
    Type Input(const std::string name=std::string(""));


    // access
    std::fstream & GetFstream() { return m_file; }
    
    unsigned long Nx() const { return m_nx; }
    unsigned long Ny() const { return m_ny; }


    // helpers
    void FillIn(const std::string &);
    void Shorten(std::string&);
  };
}

#endif //  ATOOLS

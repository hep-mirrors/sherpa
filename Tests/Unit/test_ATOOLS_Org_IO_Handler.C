#include <catch2/catch_all.hpp>
#include "ATOOLS/Org/IO_Handler.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include <cstdio>
#include <fstream>
#include <string>

using namespace ATOOLS;

namespace {
  // writes a 2x2 matrix in the layout of AMEGIC's colour matrix files and
  // reads it back with IO_Handler::MatrixInput<Complex>
  Complex** ReadMatrix(const std::string& row0, const std::string& row1) {
    if (!msg) msg = new Message();
    // in MPI builds My_File broadcasts the file from rank 0, so reading it
    // needs MPI and the mpi global, which only Sherpa's main sets up
    if (!mpi) {
#ifdef USING__MPI
      int initialised(0);
      MPI_Initialized(&initialised);
      if (!initialised) MPI_Init(NULL,NULL);
#endif
      mpi = new My_MPI();
    }
    const std::string file("test_ATOOLS_Org_IO_Handler.dat");
    {
      std::ofstream out(file);
      out<<"[2;2]{{"<<row0<<"};\n{"<<row1<<"}}\n";
    }
    IO_Handler ioh;
    REQUIRE(ioh.SetFileNameRO(file));
    Complex** m(NULL);
    try { m=ioh.MatrixInput<Complex>("",2,2); }
    catch (...) { std::remove(file.c_str()); throw; }
    std::remove(file.c_str());
    return m;
  }
}

TEST_CASE("IO_Handler MatrixInput", "[ATOOLS::Org::IO_Handler][ATOOLS][Org]") {
  SECTION("std format") {
    Complex** m(ReadMatrix("(12,-0);(0,6)", "(0,-6);(5.5,0)"));
    CHECK(m[0][0]==Complex(12.,0.));
    CHECK(m[0][1]==Complex(0.,6.));
    CHECK(m[1][0]==Complex(0.,-6.));
    CHECK(m[1][1]==Complex(5.5,0.));
  }
  SECTION("real numbers") {
    Complex** m(ReadMatrix("12;0", "-1.5;5.5"));
    CHECK(m[0][0]==Complex(12.,0.));
    CHECK(m[1][0]==Complex(-1.5,0.));
  }
  SECTION("unparsable entry") {
    REQUIRE_THROWS(ReadMatrix("12;i*6", "i*-6;5.5"));
  }
  SECTION("partially parsable entry") {
    REQUIRE_THROWS(ReadMatrix("12;0+i*6", "0-i*6;5.5"));
  }
}

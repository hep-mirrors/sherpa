#include "ATOOLS/Org/Gzip_Stream.H"
#include <iostream>
#include <string.h>  // for memcpy

#include "Data_Reader.H"
#include "Exception.H"

namespace ATOOLS {

// --------------------------------------
// class Gzip_Stream
// --------------------------------------

  Gzip_Stream::Gzip_Stream() : p_stream(NULL), m_streamtype(0) {
    Data_Reader dr(" ",";","!","=");
    int use_gzip = dr.GetValue<int>("USE_GZIP",
#ifdef USING__GZIP
                                    1
#else
                                    0
#endif
                                    );

    if (use_gzip==0) {
      p_stream = new std::ofstream();
      m_streamtype=1;
    }
    else {
#ifdef USING__GZIP
      p_stream = new ogzstream();
      m_streamtype=2;
#else
      THROW(fatal_error, "Asked for GZIP but did not compile with GZIP support.");
#endif
    }
  }

  Gzip_Stream::~Gzip_Stream() {
    delete p_stream;
  }

  std::ostream* Gzip_Stream::stream() {
    return p_stream;
  }

  void Gzip_Stream::open(const std::string &name)
  {
    if (m_streamtype==1) {
      std::ofstream* mystream = dynamic_cast<std::ofstream*>(p_stream);
      if (mystream==NULL) THROW(fatal_error, "Internal error 1");
      mystream->open(name.c_str());
    }
    if (m_streamtype==2) {
#ifdef USING__GZIP
      ATOOLS::ogzstream* mystream = dynamic_cast<ATOOLS::ogzstream*>(p_stream);
      if (mystream==NULL) THROW(fatal_error, "Internal error 2");
      mystream->open(name.c_str());
#else
      THROW(fatal_error, "Asked for GZIP but did not compile with GZIP support.");
#endif
    }
    return;
  }

  void Gzip_Stream::close()
  {
    if (m_streamtype==1) {
      std::ofstream* mystream = dynamic_cast<std::ofstream*>(p_stream);
      if (mystream==NULL) THROW(fatal_error, "Internal error 1");
      mystream->close();
    }
    if (m_streamtype==2) {
#ifdef USING__GZIP
      ATOOLS::ogzstream* mystream = dynamic_cast<ATOOLS::ogzstream*>(p_stream);
      if (mystream==NULL) THROW(fatal_error, "Internal error 2");
      mystream->close();
#else
      THROW(fatal_error, "Asked for GZIP but did not compile with GZIP support.");
#endif
    }
    return;
  }


#ifdef USING__GZIP

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See header file for user classes.
// ----------------------------------------------------------------------------

// --------------------------------------
// class gzstreambuf:
// --------------------------------------

gzstreambuf* gzstreambuf::open( const char* name, int open_mode) {
    if ( is_open())
        return (gzstreambuf*)0;
    mode = open_mode;
    // no append nor read/write mode
    if ((mode & std::ios::ate) || (mode & std::ios::app)
        || ((mode & std::ios::in) && (mode & std::ios::out)))
        return (gzstreambuf*)0;
    char  fmode[10];
    char* fmodeptr = fmode;
    if ( mode & std::ios::in)
        *fmodeptr++ = 'r';
    else if ( mode & std::ios::out)
        *fmodeptr++ = 'w';
    *fmodeptr++ = 'b';
    *fmodeptr = '\0';
    file = gzopen( name, fmode);
    if (file == 0)
        return (gzstreambuf*)0;
    opened = 1;
    return this;
}

gzstreambuf * gzstreambuf::close() {
    if ( is_open()) {
        sync();
        opened = 0;
        if ( gzclose( file) == Z_OK)
            return this;
    }
    return (gzstreambuf*)0;
}

int gzstreambuf::underflow() { // used for input buffer only
    if ( gptr() && ( gptr() < egptr()))
        return * reinterpret_cast<unsigned char *>( gptr());

    if ( ! (mode & std::ios::in) || ! opened)
        return EOF;
    // Josuttis' implementation of inbuf
    int n_putback = gptr() - eback();
    if ( n_putback > 4)
        n_putback = 4;
    memcpy( buffer + (4 - n_putback), gptr() - n_putback, n_putback);

    int num = gzread( file, buffer+4, bufferSize-4);
    if (num <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    setg( buffer + (4 - n_putback),   // beginning of putback area
          buffer + 4,                 // read position
          buffer + 4 + num);          // end of buffer

    // return next character
    return * reinterpret_cast<unsigned char *>( gptr());    
}

int gzstreambuf::flush_buffer() {
    // Separate the writing of the buffer from overflow() and
    // sync() operation.
    int w = pptr() - pbase();
    if ( gzwrite( file, pbase(), w) != w)
        return EOF;
    pbump( -w);
    return w;
}

int gzstreambuf::overflow( int c) { // used for output buffer only
    if ( ! ( mode & std::ios::out) || ! opened)
        return EOF;
    if (c != EOF) {
        *pptr() = c;
        pbump(1);
    }
    if ( flush_buffer() == EOF)
        return EOF;
    return c;
}

int gzstreambuf::sync() {
    // Changed to use flush_buffer() instead of overflow( EOF)
    // which caused improper behavior with std::endl and flush(),
    // bug reported by Vincent Ricard.
    if ( pptr() && pptr() > pbase()) {
        if ( flush_buffer() == EOF)
            return -1;
    }
    return 0;
}

// --------------------------------------
// class gzstreambase:
// --------------------------------------

gzstreambase::gzstreambase( const char* name, int mode) {
    init( &buf);
    open( name, mode);
}

gzstreambase::~gzstreambase() {
    buf.close();
}

void gzstreambase::open( const char* name, int open_mode) {
    if ( ! buf.open( name, open_mode))
        clear( rdstate() | std::ios::badbit);
}

void gzstreambase::close() {
    if ( buf.is_open())
        if ( ! buf.close())
            clear( rdstate() | std::ios::badbit);
}

#endif
  
} // namespace ATOOLS


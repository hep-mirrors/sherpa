#include "ATOOLS/Org/RUsage.H"

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif
#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

size_t ATOOLS::GetPeakRSS( )
{
#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  /* BSD, Linux, and OSX -------------------------------------- */
  struct rusage rusage;
  getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
  return (size_t)rusage.ru_maxrss;
#else
  return (size_t)(rusage.ru_maxrss * 1024L);
#endif
#else
  return (size_t)0L;          /* Unsupported. */
#endif
}

size_t ATOOLS::GetCurrentRSS( )
{
#if defined(__APPLE__) && defined(__MACH__)
  /* OSX ------------------------------------------------------ */
  struct mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		  (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
    return (size_t)0L;      /* Can't access? */
  return (size_t)info.resident_size;
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
  /* Linux ---------------------------------------------------- */
  long rss = 0L;
  FILE* fp = NULL;
  if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
    return (size_t)0L;      /* Can't open? */
  if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
      fclose( fp );
      return (size_t)0L;      /* Can't read? */
    }
  fclose( fp );
  return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);
#else
  return (size_t)0L;          /* Unsupported. */
#endif
}

#ifndef ATOOLS_Math_Table_Cache_H
#define ATOOLS_Math_Table_Cache_H

#include "ATOOLS/Org/My_File.H"

#include <fstream>
#include <memory>
#include <sstream>
#include <string>

namespace ATOOLS {

  // Versioned, MPI-safe disk cache for the look-up tables in Lookup_Tables.H
  // (OneDim_Table / TwoDim_Table) -- in fact for any type providing
  //     void Write(std::ostream&) const;
  //     static std::unique_ptr<T> Read(std::istream&);
  //
  // File layout (text, so it round-trips through My_File's rank-0-read +
  // broadcast):
  //     <marker> <format-version>\n
  //     <key>\n
  //     <table body written by Table::Write>
  //     <end-marker>\n
  // The key string fully identifies the cached content; the caller is
  // responsible for encoding every relevant parameter (and any content version)
  // into it. Load/Save go through My_File, so only rank 0 ever touches the disk
  // and no shared filesystem is assumed across ranks.
  namespace Table_Cache {

    // File-format framing constants. Bump s_format_version if this layout
    // (as opposed to the caller's content) ever changes.
    constexpr const char* s_marker = "TABLECACHE";
    constexpr const char* s_end_marker = "TABLECACHE_END";
    constexpr int s_format_version = 1;

    // Returns the cached table, or nullptr on a miss: file absent, marker /
    // format-version / key mismatch, parse error, or a torn/truncated write
    // (detected via the trailing end marker). MPI-collective: My_In_File has
    // rank 0 read the file and broadcast its content to every rank (a missing
    // file makes Open() return false on all ranks alike), so this MUST be
    // reached by all ranks together.
    template <class Table>
    std::unique_ptr<Table> Load(const std::string& path, const std::string& key)
    {
      My_In_File in(path);
      if (!in.Open()) return nullptr;
      std::istringstream is(in.FileContent());
      std::string headerline;
      if (!std::getline(is, headerline)) return nullptr;
      std::istringstream header(headerline);
      std::string marker;
      int version = 0;
      header >> marker >> version;
      if (!header || marker != s_marker || version != s_format_version)
        return nullptr;
      std::string fkey;
      if (!std::getline(is, fkey) || fkey != key) return nullptr;
      try {
        auto table = Table::Read(is);
        // Require the trailing end marker: a torn/truncated write (or a body
        // with extra trailing tokens) leaves it missing/mismatched, so we treat
        // the file as a miss rather than load a partial table.
        std::string endmarker;
        if (!(is >> endmarker) || endmarker != s_end_marker) return nullptr;
        return table;
      } catch (...) {
        return nullptr; // corrupt / truncated cache: recompute
      }
    }

    // Writes the table to disk through My_Out_File (rank 0 only; other ranks
    // build the content and discard it on close). Open() merely sets up the
    // in-memory buffer (it does not fail for an output file); the write is
    // best-effort, so a failure to flush to disk is non-fatal and simply makes
    // the next Load a miss.
    template <class Table>
    void Save(const std::string& path, const std::string& key,
              const Table& table)
    {
      My_Out_File out(path);
      out.Open();
      *out << s_marker << " " << s_format_version << "\n" << key << "\n";
      table.Write(*out);
      *out << s_end_marker << "\n";
    }

  } // namespace Table_Cache

} // namespace ATOOLS

#endif

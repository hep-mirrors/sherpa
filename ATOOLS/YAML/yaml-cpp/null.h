#ifndef ATOOLS_YAML_NULL_H_62B23520_7C8E_11DE_8A39_0800200C9A66
#define ATOOLS_YAML_NULL_H_62B23520_7C8E_11DE_8A39_0800200C9A66

#if defined(_MSC_VER) ||                                            \
    (defined(__GNUC__) && (__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || \
     (__GNUC__ >= 4))  // GCC supports "pragma once" correctly since 3.4
#pragma once
#endif

#include "ATOOLS/YAML/yaml-cpp/dll.h"
#include <string>

namespace SHERPA_YAML {
class Node;

struct YAML_CPP_API _Null {};
inline bool operator==(const _Null&, const _Null&) { return true; }
inline bool operator!=(const _Null&, const _Null&) { return false; }

YAML_CPP_API bool IsNull(const Node& node);  // old API only
YAML_CPP_API bool IsNullString(const std::string& str);

extern YAML_CPP_API _Null Null;
}

#endif  // ATOOLS_YAML_NULL_H_62B23520_7C8E_11DE_8A39_0800200C9A66

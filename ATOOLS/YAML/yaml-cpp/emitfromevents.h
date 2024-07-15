#ifndef ATOOLS_YAML_EMITFROMEVENTS_H_62B23520_7C8E_11DE_8A39_0800200C9A66
#define ATOOLS_YAML_EMITFROMEVENTS_H_62B23520_7C8E_11DE_8A39_0800200C9A66

#if defined(_MSC_VER) ||                                            \
    (defined(__GNUC__) && (__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || \
     (__GNUC__ >= 4))  // GCC supports "pragma once" correctly since 3.4
#pragma once
#endif

#include <stack>

#include "ATOOLS/YAML/yaml-cpp/anchor.h"
#include "ATOOLS/YAML/yaml-cpp/emitterstyle.h"
#include "ATOOLS/YAML/yaml-cpp/eventhandler.h"

namespace SHERPA_YAML {
struct Mark;
}  // namespace SHERPA_YAML

namespace SHERPA_YAML {
class Emitter;

class EmitFromEvents : public EventHandler {
 public:
  EmitFromEvents(Emitter& emitter);

  void OnDocumentStart(const Mark& mark) override;
  void OnDocumentEnd() override;

  void OnNull(const Mark& mark, anchor_t anchor) override;
  void OnAlias(const Mark& mark, anchor_t anchor) override;
  void OnScalar(const Mark& mark, const std::string& tag,
                        anchor_t anchor, const std::string& value) override;

  void OnSequenceStart(const Mark& mark, const std::string& tag,
                               anchor_t anchor, EmitterStyle::value style) override;
  void OnSequenceEnd() override;

  void OnMapStart(const Mark& mark, const std::string& tag,
                          anchor_t anchor, EmitterStyle::value style) override;
  void OnMapEnd() override;

 private:
  void BeginNode();
  void EmitProps(const std::string& tag, anchor_t anchor);

 private:
  Emitter& m_emitter;

  struct State {
    enum value { WaitingForSequenceEntry, WaitingForKey, WaitingForValue };
  };
  std::stack<State::value> m_stateStack;
};
}

#endif  // ATOOLS_YAML_EMITFROMEVENTS_H_62B23520_7C8E_11DE_8A39_0800200C9A66

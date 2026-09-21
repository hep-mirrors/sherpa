add_subdirectory(Amplitude)
add_subdirectory(DipoleSubtraction)
add_subdirectory(Main)
add_subdirectory(Phasespace)
add_subdirectory(String)
add_subdirectory(Amplitude/Zfunctions)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(AmegicCore SHARED
    $<TARGET_OBJECTS:Amplitude> $<TARGET_OBJECTS:DipoleSubtraction>
    $<TARGET_OBJECTS:Amegic> $<TARGET_OBJECTS:AmegicPSGen>
    $<TARGET_OBJECTS:String> $<TARGET_OBJECTS:Zfunctions>)
  target_link_libraries(AmegicCore PRIVATE Amplitude DipoleSubtraction Amegic AmegicPSGen String Zfunctions)
  sherpa_handle_shared_library(AmegicCore AmegicCore)
  add_dependencies(AmegicCore AmegicAmplitude_Git_Info DipoleSubtraction_Git_Info AmegicMain_Git_Info AmegicPhasespace_Git_Info String_Git_Info Zfunctions_Git_Info)
endif()

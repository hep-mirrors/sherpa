##### Automatically or semiautomaticaly generated, do not edit.#########
#Ladders/Makefile.am
 
set(GITTAG  SHRiMPS/Ladders)

set(Ladders_SOURCES Ladder_Particle.C
                    T_Prop.C
                    Ladder.C
                    Ladder_Generator_Base.C
                    Ladder_Generator_Eik.C
                    Ladder_Generator_QT.C
                    Ladder_Generator_KT.C
                    Ladder_Generator_Seeded.C
                    MEs.C
                    Primary_Ladders.C)

set(Ladders_HEADERS Ladder_Particle.H
                    T_Prop.H
                    Ladder.H
                    Ladder_Generator_Base.H
                    Ladder_Generator_Eik.H
                    Ladder_Generator_QT.H
                    Ladder_Generator_KT.H
                    Ladder_Generator_Seeded.H
                    MEs.H
                    Primary_Ladders.H)

sherpa_create_git_info(Ladders)
 
add_library(ShrimpsLadders SHARED ${Ladders_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsLadders ShrimpsLadders)

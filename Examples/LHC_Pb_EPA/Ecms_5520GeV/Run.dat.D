(run){
  EVENTS = 10
}(run)

(beam){
  BEAM_1          = 1000822080
  BEAM_ENERGY_1   = 2760.
  BEAM_SPECTRUM_1 = EPA

  BEAM_2          = 1000822080
  BEAM_ENERGY_2   = 2760.
  BEAM_SPECTRUM_2 = EPA  
}(beam)


(processes){
  Process : 22 22 -> 411 -411
  End process
}(processes)

(isr){
  BUNCH_1         = 22
  BUNCH_2         = 22
}(isr)


(selector){
  ET 11 10 10000
  ET -11 10 10000
}(selector)

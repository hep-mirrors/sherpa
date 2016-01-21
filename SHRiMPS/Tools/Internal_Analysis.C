  if (m_analyse) {
    m_histograms[string("N_ladder_naive")] = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_start")] = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_prim")]  = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_sec")]   = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_true")]  = new Histogram(0,0.0,25.0,25);
    m_histograms[string("B_naive")]        = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B_real")]         = new Histogram(0,0.0,25.0,50);
    m_histograms[string("N_ladder1_B")]    = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_all_B")] = new Histogram(0,0.0,25.0,25);
    m_histograms[string("B1_prim")]        = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B1_all")]         = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B2_prim")]        = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B2_all")]         = new Histogram(0,0.0,25.0,50);
  }
    if (m_analyse) {
      msg_Info()
	<<"Mean number of ladders: "
	<<"naive = "<<m_histograms[string("N_ladder_naive")]->Average()<<", "
	<<"start = "<<m_histograms[string("N_ladder_start")]->Average()<<", "
	<<"prim = "<<m_histograms[string("N_ladder_prim")]->Average()<<", "
	<<"true = "<<m_histograms[string("N_ladder_true")]->Average()<<".\n";
    }


  if (m_histograms.empty() || !m_analyse) return;
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();hit++) {
    histo = hit->second;
    name  = string("Ladder_Analysis/")+hit->first+string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histograms.clear();

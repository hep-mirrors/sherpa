  m_connectblobs = m_laddercols = m_updatecols = 0;
  msg_Info()<<"In "<<METHOD<<"(out = "<<m_output<<")\n";
  if (m_output) {
    msg_Info()<<"Errors: \n"
	      <<"   Not able to connect blobs "<<m_connectblobs<<";\n"
	      <<"   Wrong colours from ladder "<<m_laddercols<<";\n"
	      <<"   Not able to update colours in event "<<m_updatecols<<".\n";
  }

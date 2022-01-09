void run_cas_offinder(Cas_OFFinder &s, const char* chromfilename, const char* outfilename, const char* summaryfilename, int *cnum) {
	string filepath = chromfilename;
	bool issummary = strlen(summaryfilename) > 0;

	cerr << "Reading " << filepath << "..." << endl;
	
	if (read_fasta(filepath, s.m_chrnames, s.m_chrdata, s.m_chrpos)) {
		if (read_twobit(filepath, s.m_chrnames, s.m_chrdata, s.m_chrpos)) {
			cerr << "Skipping non-acceptable file " << filepath << "..." << endl;
			return;
		}
	}
	cerr << "Sending data to devices..." << endl;
	s.setChrData();
	cerr << "Chunk load started." << endl;
	while (s.loadNextChunk()) {
		// Find patterns in the chunk
		cerr << "Finding pattern in chunk #" << ++(*cnum) << "..." << endl;
		s.findPattern();
		cerr << "Comparing patterns in chunk #" << *cnum << "..." << endl;
		s.compareAll(outfilename, issummary);
		s.releaseLociinfo();
	}
	if (issummary) {
		cerr << "Writing summary file..." << endl;
		s.writeSummaryTable(summaryfilename);
	}
}
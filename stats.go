package main

type ContigStat struct {
	//name                string
	referenceLength     int
	referenecClean      int
	referenceDuplicated int
	allCalled           int
	allPassedCoverage   int
	allPassedProportion int
	allPassedConsensus  int
	qualityBreadth      int
	anySnps             int
	bestSnps            int
}

type SampleStat struct {
	//name                   string
	//aligner                string
	//snpcaller              string
	wasCalled              int
	passedCoverageFilter   int
	passedProportionFilter int
	qualityBreadth         int
	calledReference        int
	calledSnp              int
	calledDegen            int
}

type ContigStats []ContigStat

func (c ContigStats) Aggregate(ch chan []ContigStats) {
}

type SampleStats []SampleStat

func (s SampleStats) Aggregate(ch chan SampleStats /*, filepath string*/) {
	for stats := range ch {
		for i := range s {
			s[i].wasCalled += stats[i].wasCalled
			s[i].passedCoverageFilter += stats[i].passedCoverageFilter
			s[i].passedProportionFilter += stats[i].passedProportionFilter
			s[i].qualityBreadth += stats[i].qualityBreadth
			s[i].calledReference += stats[i].calledReference
			s[i].calledSnp += stats[i].calledSnp
			s[i].calledDegen += stats[i].calledDegen
		}
	}

	/*
		file, err := os.Create("sample_stats.tsv")
		if err != nil {
			log.Fatal(err)
		}

		br := bufio.NewWriter(file)
	*/
}

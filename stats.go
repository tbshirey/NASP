package main

import (
	"bufio"
	"log"
	"os"
	"strconv"
)

type Contigs struct {
	Reference  []byte
	Duplicates []byte
	Analyses   [][]byte
}

type ContigStat struct {
	name                string
	minCoverage         int
	minProportion       float64
	referenceLength     int
	referenceClean      int
	referenceDuplicated int
	allCalled           int
	allPassedCoverage   int
	allPassedProportion int
	allPassedConsensus  int
	qualityBreadth      int
	anySnps             int
	bestSnps            int
}

func NewContigStat(name string, minCoverage int, minProportion float64) *ContigStat {
	return &ContigStat{
		name:          name,
		minCoverage:   minCoverage,
		minProportion: minProportion,
	}
}

func (c ContigStat) Compare(contigs chan Contigs, sampleStats chan SampleStats, numSamples int) {
	var isDup, refCalled bool
	stats := make(SampleStats, numSamples)
	for contig := range contigs {
		//fmt.Println(contig.Analyses)
		c.referenceLength += len(contig.Reference)
		for i, refCall := range contig.Reference {
			switch refCall {
			case 'G', 'A', 'T', 'C':
				refCalled = true
			}
			isDup = contig.Duplicates != nil && contig.Duplicates[i] == '1'
			if isDup {
				c.referenceDuplicated++
			}
			for j, analysis := range contig.Analyses {
				if analysis == nil || len(analysis) < i+1 || !refCalled {
					continue
				}
				switch analysis[i] {
				default:
					stats[j].calledDegen++
				case 'X', 'N':
					break
				case refCall:
					stats[j].calledSnp++
				case 'G', 'A', 'T', 'C':
					// TODO: check pass cov/prop
					if analysis[i] == refCall {
						stats[j].calledReference++
					}

				}
			}
		}
	}
	/*
		var isDup bool
		//stats := make(SampleStats, numSamples)
		statsCh := make(chan SampleStats, 4096)
		for contig := range contigs {
			c.referenceLength += len(contig.Reference)

			for i, refCall := range contig.Reference {
				isDup = contig.Duplicates != nil && contig.Duplicates[i] == '1'
				go c.compare(statsCh, i, refCall, isDup, contig.Analyses)
			}
			//fmt.Printf("%s", ch)
		}
	*/
	/*
		ok := true
		for ok {
			_, ok = <-contigs
		}
	*/
}

func (c ContigStat) compare(statsCh chan SampleStats, i int, refCall byte, isDup bool, analyses [][]byte) {
	stats := make(SampleStats, len(analyses))
	for j, analysis := range analyses {
		if len(analysis) < i+1 {
			continue
		}
		if analysis[i] == refCall {
			stats[j].calledReference++
		}
	}
	statsCh <- stats
}

type SampleStat struct {
	//name                   string
	//aligner                string
	//snpcaller              string
	wasCalled              int64
	passedCoverageFilter   int64
	passedProportionFilter int64
	qualityBreadth         int64
	calledReference        int64
	calledSnp              int64
	calledDegen            int64
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
		statsPool.Put(stats)
	}
}

func (s SampleStats) Reset() {
	for i := range s {
		s[i].wasCalled = 0
		s[i].passedCoverageFilter = 0
		s[i].passedProportionFilter = 0
		s[i].qualityBreadth = 0
		s[i].calledReference = 0
		s[i].calledSnp = 0
		s[i].calledDegen = 0
	}
}

func (s SampleStats) WriteStats() {
	file, err := os.Create("sample_stats.tsv")
	if err != nil {
		log.Fatal(err)
	}
	bw := bufio.NewWriter(file)
	bw.Write([]byte("Sample\tSample::Analysis\twas_called\twas_called (%)\tpassed_coverage_filter\tpassed_coverage_filter (%)\tpassed_proportion_filter\tpassed_proportion_filter (%)\tquality_breadth\tquality_breadth (%)\tcalled_reference\tcalled_reference (%)\tcalled_snp\tcalled_snp (%)\tcalled_degen\tcalled_degen (%)\n\n"))

	buf := make([]byte, 100)
	for _, sample := range s {
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.wasCalled, 10))
		bw.WriteByte('\t')
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.passedCoverageFilter, 10))
		bw.WriteByte('\t')
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.passedProportionFilter, 10))
		bw.WriteByte('\t')
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.qualityBreadth, 10))
		bw.WriteByte('\t')
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.calledReference, 10))
		bw.WriteByte('\t')
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.calledSnp, 10))
		bw.WriteByte('\t')
		buf = buf[:0]
		bw.Write(strconv.AppendInt(buf, sample.calledDegen, 10))
		bw.WriteByte('\n')
	}
}

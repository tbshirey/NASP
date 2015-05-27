package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/davecheney/profile"
)

//var c uint64
var c_contig uint64
var refPath = flag.String("reference", "", "Path to the reference fasta against which samples are compared")
var dupPath = flag.String("duplicates", "", "Path to the duplicates file marking possible duplicated positions")
var minCoverage = flag.Int("coverage", 0, "Filter positions below this coverage/depth threshold")
var minProportion = flag.Float64("proportion", 0.0, "Filter positions below this proportion threshold")

var NUM_SAMPLES int

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())

	defer profile.Start(&profile.Config{
		//	CPUProfile:   true,
		MemProfile: true,
		//BlockProfile: true,
		ProfilePath: ".",
		//NoShutdownHook: true,
	}).Stop()

	t0 := time.Now()
	defer func() {
		fmt.Println(time.Now().Sub(t0))
	}()

	flag.Parse()

	if *refPath == "" {
		flag.Usage()
		os.Exit(1)
	}

	if len(flag.Args()) < 1 {
		log.Println("Warning: Expected the command to include a list or glob pattern of sample files.")
	}

	analyses := NewSampleAnalyses(flag.Args()...)

	NUM_SAMPLES = len(analyses)

	log.Println("Indexing all samples", time.Now().Sub(t0))

	reference, err := NewReference(*refPath, *dupPath)
	if err != nil {
		log.Fatal(err)
	}

	ch := make(chan chan *Position, 100)
	defer func() {
		wg.Add(1)
		close(ch)
		wg.Wait()
	}()
	go func(c chan chan *Position) {
		defer wg.Done()
		/*
			for position := range c {
				fmt.Printf("%s\n\n", position.callStr)
			}
		*/
		writeMaster(c, NUM_SAMPLES)
	}(ch)

	//	var err error
	var name string
	var ref, dup []byte
	var isPrefix bool
	pos := make([][]byte, NUM_SAMPLES)
	for {
		name, err = reference.NextContig()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		if err = analyses.SeekContig(name); err != nil {
			log.Fatal(err)
		}

		log.Println("Scanning", name, time.Now().Sub(t0))

		isPrefix = true
		for isPrefix {
			ref, dup, isPrefix, err = reference.ReadPositions()
			if err == io.EOF {
				break
			} else if err != nil {
				log.Fatalf("Error scanning reference contig %s: %s", name, err.Error())
			}
			if err = analyses.ReadPositions(pos); err != nil {
				log.Fatalf("Error reading contig %s: %s", name, err.Error())
			}
			if isPrefix {
				log.Println("Scanning...", name, time.Now().Sub(t0))
			}
			fmt.Printf("len(ref) = %d\n", len(ref))
			wg.Add(1)
			c := make(chan *Position, 100)
			ch <- c
			go analyzePositions(c, ref, dup, pos)
		}
		//fmt.Println("NumGoroutine", runtime.NumGoroutine())
	}
	wg.Wait()
}

var wg sync.WaitGroup

// Assumes all positions are uppercase
func analyzePositions(ch chan *Position, ref, dup []byte, analyses [][]byte) {
	defer func() {
		fmt.Println("Shutdown NumGoroutine", runtime.NumGoroutine())
		wg.Done()
		close(ch)
	}()
	var call byte
	stats := statsPool.Get().(SampleStats)
	for i, refCall := range ref {
		position := positionPool.Get().(*Position)

		switch refCall {
		default:
			position.pattern[0] = 'N'
		case 'G', 'A', 'T', 'C':
			position.isReferenceClean = true
			position.pattern[0] = '1'
		}
		position.callStr[0] = refCall

		position.isReferenceDuplicated = dup != nil && dup[i] == '1'

		for j := range analyses {
			if analyses[j] == nil || len(analyses[j]) < i+1 {
				//position.isAllCalled = false
				//position.callStr[(j+1)<<1-1] = '\t'
				//position.callStr[(j+1)<<1] = 'X'
				//continue
				call = 'X'
			} else {
				defer func(j, i, l int, analysis []byte) {
					if recover() != nil {
						/*
							fmt.Printf("QUOX: %s\n", analysis)
							fmt.Println("BAZQUOX", analysis)
							fmt.Println("FOOBAR:", j, i, l)
							fmt.Println("QUOXX", analyses[j][i])
							os.Exit(1)
						*/
					}
				}(j, i, len(analyses[j]), analyses[j])
				//fmt.Println(j, i, len(analyses[j]))
				call = analyses[j][i]
				//call = 'Z'
			}
			position.callStr[2*(j+1)-1] = '\t'
			position.callStr[2*(j+1)] = call

			// TODO _is_pass_filter
			stats[j].passedCoverageFilter++
			stats[j].passedProportionFilter++
			position.passedCoverage++
			position.passedProportion++
			position.passedDepthFilter[j+1] = '-'
			position.passedProportionFilter[j+1] = '-'

			// TODO: Sample consensus

			// It can be called anything so long as it was called something
			// X and N
			isDegen := false
			wasCalled := true

			switch call {
			case 'G':
				position.g++
			case 'A':
				position.a++
			case 'T':
				position.t++
			case 'C':
				position.c++
			case 'X', 'N':
				wasCalled = false
				fallthrough
			default:
				isDegen = true
				position.isAllPassConsensus = false
				position.isAllQualityBreadth = false
				position.n++
			}

			if wasCalled {
				stats[j].wasCalled++
				position.wasCalled++
				position.callWasMade[j] = 'Y'
			} else {
				position.callWasMade[j] = 'N'
				position.isAllCalled = false
				position.isAllQualityBreadth = false
			}

			// TODO: is pass cov/prop
			if wasCalled && position.isReferenceClean {
				if isDegen {
					position.calledDegen++
					if !position.isReferenceDuplicated {
						stats[j].calledDegen++
					}
				} else if refCall == call {
					position.calledReference++
					if !position.isReferenceDuplicated {
						stats[j].calledReference++
					}
				} else if !position.isReferenceDuplicated {
					// Called A/C/G/T and doesn't match the reference
					position.calledSnp++
					position.isMissingMatrix = true
				}

			} else {
				position.isAllQualityBreadth = false
			}
		}
		ch <- position
	}
	statsPool.Put(stats)
}

type Position struct {
	// General Stats
	isAllCalled           bool
	isReferenceClean      bool
	isReferenceDuplicated bool
	isAllPassCoverage     bool
	isAllPassProportion   bool
	isAllPassConsensus    bool
	isAllQualityBreadth   bool
	isBestSnp             bool
	isMissingMatrix       bool

	//all_sample_stats=all_sample_stats,

	// Missing Data Matrix condition - at least one SampleAnalysis passes quality_breadth and is a SNP.
	//is_missing_matrix=is_missing_matrix,

	// NASP Master Matrix
	// Counters
	wasCalled        int
	calledReference  int
	calledSnp        int
	calledDegen      int
	passedCoverage   int
	passedProportion int
	a                int
	c                int
	g                int
	t                int
	n                int

	// Strings
	callStr []byte
	// TODO
	//maskedCallStr          []byte
	callWasMade            []byte
	passedDepthFilter      []byte
	passedProportionFilter []byte
	pattern                []byte
}

var positionPool = sync.Pool{
	New: func() interface{} {
		// +1 for reference
		numSamples := NUM_SAMPLES + 1
		//position := &Position{
		return &Position{
			isAllCalled: true,
			//isReferenceClean:      true,
			//isReferenceDuplicated: true,
			isAllPassCoverage:   true,
			isAllPassProportion: true,
			isAllPassConsensus:  true,
			isAllQualityBreadth: true,
			isBestSnp:           true,
			//isMissingMatrix: true,

			// 2*(NUM_SAMPLES + reference) - 1
			// allows space for a \t between each call
			callStr: make([]byte, 2*(NUM_SAMPLES)+1),
			// TODO
			//maskedCallStr:          make([]byte, numSamples),
			callWasMade:            make([]byte, numSamples),
			passedDepthFilter:      make([]byte, numSamples),
			passedProportionFilter: make([]byte, numSamples),
			pattern:                make([]byte, numSamples),
		}
		/*
			for i := 1; i < len(position.callStr)-1; i += 2 {
				position.callStr[i] = '\t'
			}
			return position
		*/
	},
}

var statsPool = sync.Pool{
	New: func() interface{} {
		return make(SampleStats, NUM_SAMPLES)
	},
}

type SampleAnalysis interface {
	//Contig(done chan struct{}, name string) chan byte
	//Name() string
	//ScanPositions(name string) ([]byte, error)
	SeekContig(name string) error
	ReadPositions() ([]byte, error)
}

// SampleAnalyses implements sort.Interface for []SampleAnalysis based on
// the identifier field.
type SampleAnalyses []SampleAnalysis

func (s SampleAnalyses) ReadPositions(positions [][]byte) error {
	for i, analysis := range s {
		pos, err := analysis.ReadPositions()
		// FIXME: EOF will be common, but silently ignoring it might not be good either.
		switch err {
		default:
			return err
		case bufio.ErrBufferFull, io.EOF:
			break
		}
		positions[i] = pos
	}
	return nil
}

/**
 * SeekContig moves the io.Reader for each analysis onto the first position
 * of the contig with the given name, or to io.EOF if it does not exist.
 */
func (s SampleAnalyses) SeekContig(name string) (err error) {
	for _, analysis := range s {
		if err = analysis.SeekContig(name); err != nil {
			return err
		}
	}
	return nil
}

func (s SampleAnalyses) Len() int {
	return len(s)
}

/*
func (s SampleAnalyses) Less(i, j int) {
	return s[i].Identifier < s[j].Identifier
}
*/

func (s SampleAnalyses) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

func NewSampleAnalyses(filepaths ...string) SampleAnalyses {
	analyses := make(SampleAnalyses, len(filepaths))

	ch := make(chan SampleAnalysis, len(flag.Args()))
	defer close(ch)
	for _, path := range filepaths {
		go func(ch chan SampleAnalysis, path string) {
			switch {
			default:
				log.Fatal("Unknown sample analysis: " + path)
			case strings.HasSuffix(path, "fasta"):
				// NOTE: This will accept any path matching *fasta$
				// ex: .fasta, .frankenfasta, notReallyAfasta
				fasta, err := NewFasta(path, true)
				if err != nil {
					log.Fatal(err)
				}
				ch <- fasta
				//case strings.HasSuffix(path, "vcf"):
				//	ch <- NewVcf(path)
			}
			return
		}(ch, path)
	}

	for i := range analyses {
		analyses[i] = <-ch
	}

	return analyses
}

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

const (
	// defaultBufSize is the maximum number of calls read from a contig at a time.
	// defaultBufSize * numFiles * maxQueuedBlocks roughly translates to RAM usage.
	defaultBufSize  = 4096
	maxQueuedBlocks = 10
)

var NUM_SAMPLES int
var t0 = time.Now()

var refPath = flag.String("reference", "", "Path to the reference fasta against which samples are compared")
var dupPath = flag.String("duplicates", "", "Path to the duplicates file marking possible duplicated positions")
var minCoverage = flag.Int("coverage", 0, "Filter positions below this coverage/depth threshold")
var minProportion = flag.Float64("proportion", 0.0, "Filter positions below this proportion threshold")

func main() {
	var wg sync.WaitGroup
	runtime.GOMAXPROCS(runtime.NumCPU())

	// Begin Development Profiling
	defer profile.Start(&profile.Config{
		//	CPUProfile:   true,
		MemProfile:   true,
		BlockProfile: true,
		ProfilePath:  ".",
		//NoShutdownHook: true,
	}).Stop()

	t0 := time.Now()
	defer func() {
		fmt.Println(time.Now().Sub(t0))
	}()
	// End Development Profiling

	flag.Parse()

	if *refPath == "" {
		flag.Usage()
		os.Exit(1)
	}

	if len(flag.Args()) < 1 {
		log.Println("Warning: Expected the command to include a list or glob pattern of sample files.")
	}

	analyses := NewSampleAnalyses(flag.Args()...)
	log.Println("Samples indexed", time.Now().Sub(t0))

	NUM_SAMPLES = len(analyses)

	reference, err := NewReference(*refPath, *dupPath)
	if err != nil {
		log.Fatal(err)
	}

	// Queue limits the rate positions are read from the files and ensures the
	// chunks are assembled in the correct order.
	queue := make(chan Chunk, maxQueuedBlocks)
	statsChan := make(chan SampleStats, maxQueuedBlocks)
	defer func() {
		wg.Add(1)
		close(queue)
		wg.Wait()
		wg.Add(1)
		close(statsChan)
		wg.Wait()
	}()

	go func(queue chan Chunk) {
		defer wg.Done()
		// queue must be closed to exit
		writeMaster(queue, NUM_SAMPLES)
	}(queue)

	go func(statsChan chan SampleStats) {
		defer wg.Done()
		// statsCh must be closed to exit
		stats := make(SampleStats, NUM_SAMPLES)
		stats.Aggregate(statsChan)
		stats.WriteStats()
	}(statsChan)

	readPositions(queue, statsChan, reference, analyses)
}

type Chunk struct {
	contigName    string
	positionsChan chan []*Position
}

// Read positions from the reference and sample files
func readPositions(queue chan Chunk, statsChan chan SampleStats, reference *Reference, analyses SampleAnalyses) {
	var err error
	var name string
	var ref, dup []byte
	var calls [][]byte
	var isPrefix bool

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
			ref, dup, isPrefix, err = reference.ReadPositions(defaultBufSize)
			if err == io.EOF {
				break
			} else if err != nil {
				log.Fatalf("Error scanning reference contig %s: %s", name, err.Error())
			}

			// TODO: cov / prop [][]bool
			// calls is a chunk of calls that can be min(remaining sample contig length, defaultBufSize)
			calls, err = analyses.ReadPositions(defaultBufSize)
			if err != nil {
				log.Fatalf("Error reading contig %s: %s", name, err.Error())
			}
			if isPrefix {
				log.Println("Scanning...", name, time.Now().Sub(t0))
			}

			chunk := Chunk{
				contigName:    name,
				positionsChan: make(chan []*Position, 1),
			}

			queue <- chunk

			go analyzePositions(chunk.positionsChan, statsChan, ref, dup, calls)
		}
	}
}

// TODO: rename ch to positionsChan
func analyzePositions(ch chan []*Position, statsChan chan SampleStats, ref, dup []byte, analyses [][]byte) {
	defer func() {
		fmt.Println("Shutdown NumGoroutine", runtime.NumGoroutine())
		close(ch)
	}()

	var call byte
	stats := statsPool.Get().(SampleStats)
	positions := positionsPool.Get().([]*Position)[:len(ref)]

	for i := range ref {
		refCall := ToUpper(ref[i])
		position := positionPool.Get().(*Position)

		switch refCall {
		default:
			position.pattern[0] = 'N'
		case 'G', 'A', 'T', 'C':
			position.isReferenceClean = true
			position.pattern[0] = '1'
		}
		position.callStr[0] = refCall

		position.isReferenceDuplicated = len(dup) > 0 && dup[i] == '1'

		for j := range analyses {
			if analyses[j] == nil || len(analyses[j]) < i+1 {
				call = 'X'
			} else {
				call = ToUpper(analyses[j][i])
			}
			//position.callStr[2*(j+1)-1] = '\t'
			//position.callStr[2*(j+1)] = call
			position.callStr[j+1] = call

			// TODO _is_pass_filter
			stats[j].passedCoverageFilter++
			stats[j].passedProportionFilter++
			position.passedCoverage++
			position.passedProportion++
			position.passedDepthFilter[j] = '-'
			position.passedProportionFilter[j] = '-'

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
				position.n++
				isDegen = true
				position.isAllPassConsensus = false
				position.isAllQualityBreadth = false
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
		positions[i] = position
	}
	callsPool.Put(analyses)
	ch <- positions
	statsChan <- stats
}

type SampleAnalysis interface {
	//Name() string
	SeekContig(name string) error
	ReadPositions(n int) ([]byte, error)
}

// SampleAnalyses implements sort.Interface for []SampleAnalysis based on
// the identifier field.
type SampleAnalyses []SampleAnalysis

func (s SampleAnalyses) ReadPositions(n int) ([][]byte, error) {
	positions := callsPool.Get().([][]byte)
	for i, analysis := range s {
		pos, err := analysis.ReadPositions(n)
		// FIXME: EOF will be common, but silently ignoring it might not be good either.
		switch err {
		default:
			positionsPool.Put(positions)
			return nil, err
		case bufio.ErrBufferFull, io.EOF:
			break
		case nil:
			positions[i] = positions[i][:len(pos)]
			copy(positions[i], pos)
		}
	}
	return positions, nil
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
				log.Fatalf("Unknown sample analysis: %s\n", path)
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

/**
 * ToUpper is reduced from the standard library version which includes unicode
 * support and consequently additional memory allocation.
 * It promotes an ASCII byte to uppercase.
 */
func ToUpper(b byte) byte {
	if 'a' <= b && b <= 'z' {
		b -= 'a' - 'A'
	}
	return b
}

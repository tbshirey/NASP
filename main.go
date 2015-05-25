package main

import (
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

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())

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

	flag.Parse()

	if *refPath == "" {
		flag.Usage()
		os.Exit(1)
	}

	if len(flag.Args()) < 1 {
		log.Println("Warning: Expected the command to include a list or glob pattern of sample files.")
	}

	analyses := NewSampleAnalyses(flag.Args()...)

	log.Println("Indexing all samples", time.Now().Sub(t0))

	reference, err := NewReference(*refPath, *dupPath)
	if err != nil {
		log.Fatal(err)
	}

	var wg sync.WaitGroup

	wg.Add(1)
	sampleStatsCh := make(chan SampleStats)
	sampleStats := make(SampleStats, len(analyses))
	go func(ch chan SampleStats) {
		defer wg.Done()
		// Will write stats to disk and shutdown when the channel is closed.
		sampleStats.Aggregate(ch)
	}(sampleStatsCh)

	var contig *ReferenceContig
	var ref, dup []byte
	var isPrefix bool

	for {
		contig, err = reference.NextContig()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		go contig.Compare(analyses, sampleStatsCh)

		isPrefix = true
		for err == nil || isPrefix {
			ref, dup, isPrefix, err = reference.PositionSlice()
			contig.RefChannel <- ref
			contig.DupChannel <- dup
		}
		//contigStatsCh <- contig
		fmt.Println(contig)
	}
	close(sampleStatsCh)

	/*
		name, err := reference.NextContig()
		if err != nil {
			log.Fatal(err)
		}

		for _, analysis := range analyses {
			line, err := analysis.ScanPositions(name)
			if err != nil {
				log.Fatal(err)
			}
			fmt.Printf("%s: %s\n", analysis.Name(), line)
		}
	*/

	// TODO: Sort analyses by Identifier
	/*
		analyses := NewSampleAnalyses(flag.Args()...)

		fmt.Println("Indexing all samples", time.Now().Sub(t0))

		positions := make([]chan byte, len(flag.Args()))

		file, err := os.Open(*refPath)
		if err != nil {
			log.Fatal(err)
		}
		br := bufio.NewReader(file)
	*/

	/*
		ch := make(chan []byte, 100)
		d := make(chan struct{})
		defer func() {
			log.Println("Close call string channel")
			close(ch)
			log.Println("Wait for all write goroutines to shutdown")
			<-d
			log.Println("Write goroutines acknowledged shutdown")
		}()
		go writeMaster(d, ch, len(positions))
	*/

	/*
		var b byte
		var r rune
		var name string
		done := make(chan struct{})
		row := make([]byte, len(positions)+1)
		var l int
		for {
			b, err = br.ReadByte()
			if err == io.EOF {
				close(done)
				fmt.Println(" Length:", l)
				return
			}
			if err != nil {
				log.Fatal(err)
			}
			r = rune(b)
			switch {
			default:
				log.Fatalf("Unexpected character in reference fasta: %c\n", b)
			case unicode.IsLetter(r):
				l++
				r = unicode.ToUpper(r)
				row[0] = byte(r)
				for i, position := range positions {
					if p, ok := <-position; ok {
						row[i+1] = p
					} else {
						row[i+1] = 'X'
					}
				}
				ch <- row
				//fmt.Printf("%s\n", row)
			case unicode.IsSpace(r):
				continue
			case r == '>':
				line, err := br.ReadSlice('\n')
				if err != nil {
					log.Fatal(err)
				}

				//line = bytes.TrimPrefix(line, []byte("franken::"))

				if idx := bytes.IndexAny(line, " \n"); idx < 0 {
					name = string(line)
				} else {
					name = string(line[:idx])
				}

				fmt.Println(" Length:", l)
				l = 0
				c_contig++
				t := time.Now().Sub(t0)
				fmt.Print(c_contig, " Scanning ", name, " ", t)

				close(done)
				done = make(chan struct{})
				for i, analysis := range analyses {
					positions[i] = analysis.Contig(done, name)
				}
			}
		}
	*/
}

type SampleAnalysis interface {
	Contig(done chan struct{}, name string) chan byte
	Name() string
	ScanPositions(name string) ([]byte, error)
}

// SampleAnalyses implements sort.Interface for []SampleAnalysis based on
// the identifier field.
type SampleAnalyses []SampleAnalysis

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

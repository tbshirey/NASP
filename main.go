package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/signal"
	"runtime"
	"strings"
	"time"
	"unicode"

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

	// TODO: Sort analyses by Identifier
	analyses := NewSampleAnalyses(flag.Args()...)

	fmt.Println("Indexing all samples", time.Now().Sub(t0))

	positions := make([]chan byte, len(flag.Args()))

	file, err := os.Open(*refPath)
	if err != nil {
		log.Fatal(err)
	}
	br := bufio.NewReader(file)

	ch := make(chan []byte, 100)
	d := make(chan struct{})
	defer func() {
		log.Println("Close call string channel")
		close(ch)
		log.Println("Wait for all write goroutines to shutdown")
		<-d
		log.Println("Write goroutines acknowledged shutdown")
	}()
	go func(done chan struct{}, ch chan []byte, numSamples int) {
		/*
					'LocusID': "{0}::{1}".format(contig_name, position),
			        #         'Reference': row.call_str[0],
			        #         '#SNPcall': row.called_snp,
			        #         # TODO: replace with n/a
			        #         '#Indelcall': '0',
			        #         '#Refcall': row.called_reference,
			        #         '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
			        #         '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
			        #         '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
			        #         '#A': row.num_A,
			        #         '#C': row.num_C,
			        #         '#G': row.num_G,
			        #         '#T': row.num_T,
			        #         # TODO: replace with n/a
			        #         '#Indel': '0',
			        #         '#NXdegen': row.num_N,
			        #         'Contig': contig_name,
			        #         'Position': position,
			        #         'InDupRegion': row.is_reference_duplicated,
			        #         'SampleConsensus': row.is_all_passed_consensus,
			        #         'CallWasMade': row.CallWasMade,
			        #         'PassedDepthFilter': row.PassedDepthFilter,
			        #         'PassedProportionFilter': row.PassedProportionFilter,
			        #         'Pattern': "".join(row.Pattern)
		*/
		c := make(chan os.Signal, 1)
		signal.Notify(c, os.Interrupt)

		var buf bytes.Buffer
		file, err := os.Create("master.tsv.gz")
		if err != nil {
			log.Fatal(err)
		}
		// Compressed Writer
		cw := gzip.NewWriter(file)
		defer func() {
			log.Println("master.tsv.gz goroutine shutting down")
			if err := cw.Close(); err != nil {
				log.Fatal(err)
			}
			if err = file.Close(); err != nil {
				log.Fatal(err)
			}
			done <- struct{}{}
		}()
		if _, err := cw.Write([]byte("LocusID\tReference\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\n")); err != nil {
			log.Fatal(err)
		}
		for i := 1; ; i++ {
			select {
			case callStr, ok := <-ch:
				if !ok {
					return
				}
				buf.Reset()
				buf.Write(callStr)
				if _, err = buf.WriteTo(cw); err != nil {
					log.Fatal(err)
				}
				if _, err = cw.Write([]byte("\n")); err != nil {
					log.Fatal(err)
				}
			case <-c:
				return
			}
		}
	}(d, ch, len(positions))

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
			fmt.Print(c_contig /*c_contig/(uint64(t/time.Second)),*/, " Scanning ", name, " ", t)

			close(done)
			done = make(chan struct{})
			for i, analysis := range analyses {
				positions[i] = analysis.Contig(done, name)
			}
		}
	}
}

type SampleAnalysis interface {
	Contig(done chan struct{}, name string) chan byte
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

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
		close(ch)
		<-d
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
		/*
			var g, a, t, c, n int
			row := make([][]byte, numSamples+22)
		*/
		c := make(chan os.Signal, 1)
		signal.Notify(c, os.Interrupt)

		//var buf bytes.Buffer
		//callStr := make([]byte, numSamples)
		//var ok bool
		file, err := os.Create("master.tsv.gz")
		if err != nil {
			log.Fatal(err)
		}
		cw := gzip.NewWriter(file)
		defer func() {
			log.Println("Write goroutine shutting down")
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
				if _, err = cw.Write([]byte(string(callStr))); err != nil {
					log.Fatal(err)
				}
				/*
					if _, err = cw.Write([]byte("LocusID\tReference\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\n")); err != nil {
						log.Fatal(err)
					}
				*/
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
	defer close(done)
	row := make([]byte, len(positions)+1)
	var l int
	for {
		b, err = br.ReadByte()
		if err == io.EOF {
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
				ch <- NewFasta(path)
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

type SampleAnalysis interface {
	Contig(done chan struct{}, name string) chan byte
}

/*
type Position struct {
	Call       byte
	Coverage   int
	Proportion float64
}
*/

type Fasta struct {
	path    string
	contigs map[string]uint
	rd      *os.File
	br      *bufio.Reader
}

func NewFasta(path string) *Fasta {
	file, err := os.Open(path)
	if err != nil {
		log.Fatal(err)
	}

	fasta := &Fasta{
		path:    path,
		contigs: make(map[string]uint),
		rd:      file,
		br:      bufio.NewReader(file),
	}

	fasta.index()
	return fasta
}

/*
func (f Fasta) Contigs() []string {
	keys := make([]string, len(f.contigs))
	i := 0
	for name := range f.contigs {
		keys[i] = name
		i += 1
	}

	return keys
}
*/

func (f Fasta) Contig(done chan struct{}, name string) chan byte {
	ch := make(chan byte, 500)
	go func(done chan struct{}, ch chan byte) {
		defer close(ch)
		var err error
		var b byte
		var r rune
		filePosition, ok := f.contigs[name]
		if !ok {
			return
		}

		_, err = f.rd.Seek(int64(filePosition), os.SEEK_SET)
		if err != nil {
			log.Fatal(err)
		}

		f.br.Reset(f.rd)
		for err != io.EOF {
			b, err = f.br.ReadByte()
			if err != nil {
				log.Fatal(err)
			}
			r = rune(b)
			switch {
			case unicode.IsSpace(r):
				continue
			case r == '>':
				return
			}

			select {
			case <-done:
				return
			case ch <- byte(unicode.ToUpper(r)):
			}
		}
	}(done, ch)

	return ch
}

/**
 * index maps the starting file position of each contig in the file
 */
func (f Fasta) index() {
	var line []byte
	var err error
	var position uint

	for err != io.EOF {
		line, err = f.br.ReadSlice('>')
		switch err {
		default:
			log.Fatal(err)
		case bufio.ErrBufferFull, io.EOF:
			break
		case nil:
			position += uint(len(line))
			line, err = f.br.ReadBytes('\n')
			if err != nil {
				log.Fatal(err)
			}
			position += uint(len(line))

			line = bytes.TrimPrefix(line, []byte("franken::"))

			var name string
			if idx := bytes.IndexAny(line, " \n"); idx < 0 {
				name = string(line)
			} else {
				name = string(line[:idx])
			}
			f.contigs[name] = position
		}
	}
}

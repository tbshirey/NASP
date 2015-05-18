package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"
	"unicode"
)

//var c uint64
var c_contig uint64
var refPath = flag.String("reference", "", "Path to the reference fasta against which samples are compared")
var dupPath = flag.String("duplicates", "", "Path to the duplicates file marking possible duplicated positions")

var EMPTY_FASTA_POSITION = &Position{
	Call:       'X',
	Coverage:   -1,
	Proportion: -1.0,
}

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())

	cprof, err := os.Create("cpu.profile")
	if err != nil {
		log.Fatal(err)
	}
	defer cprof.Close()
	pprof.StartCPUProfile(cprof)
	defer pprof.StopCPUProfile()

	mprof, err := os.Create("mem.profile")
	if err != nil {
		log.Fatal(err)
	}
	defer mprof.Close()
	defer pprof.WriteHeapProfile(mprof)

	flag.Parse()

	if *refPath == "" {
		flag.Usage()
		os.Exit(1)
	}

	if len(flag.Args()) < 1 {
		log.Println("Warning: Expected the command to include a list or glob pattern of sample files.")
	}

	t0 := time.Now()

	//reference := NewSampleAnalyses(*refPath, *dupPath)
	// TODO: Sort analyses by Identifier
	analyses := NewSampleAnalyses(flag.Args()...)
	//fmt.Println(analyses)

	fmt.Println("NewSampleAnalyses", time.Now().Sub(t0))

	positions := make([]chan *Position, len(flag.Args()))

	file, err := os.Open(*refPath)
	if err != nil {
		log.Fatal(err)
	}
	br := bufio.NewReader(file)

	var b byte
	var r rune
	var name string
	done := make(chan struct{})
	var l int
	for {
		b, err = br.ReadByte()
		if err == io.EOF {
			close(done)
			for _, p := range positions {
				<-p
				//close(p)
			}
			fmt.Println(" Length:", l)
			os.Exit(0)
		}
		if err != nil {
			log.Fatal(err)
		}
		r = unicode.ToUpper(rune(b))
		switch {
		default:
			log.Fatalf("Unexpected character in reference fasta: %c\n", r)
		case unicode.IsLetter(r):
			l++
			for _, position := range positions {
				<-position
				//p := <-position
				//fmt.Printf("%c", p.Call)
			}
			//fmt.Println()
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
			fmt.Print(c_contig, " Scanning ", name, " ", time.Now().Sub(t0))

			close(done)
			done = make(chan struct{})
			for i, analysis := range analyses {
				if positions[i] != nil {
					<-positions[i]
				}
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
	Contig(done chan struct{}, name string) chan *Position
}

type Position struct {
	Call       rune
	Coverage   int
	Proportion float64
}

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

func (f Fasta) Contigs() []string {
	keys := make([]string, len(f.contigs))
	i := 0
	for name := range f.contigs {
		keys[i] = name
		i += 1
	}

	return keys
}

func (f Fasta) emptyContig(done chan struct{}, ch chan *Position) {
	for {
		select {
		case <-done:
			return
		default:
			ch <- EMPTY_FASTA_POSITION
		}
	}
}

func (f Fasta) Contig(done chan struct{}, name string) chan *Position {
	ch := make(chan *Position, 500)
	go func(ch chan *Position) {
		var err error
		var b byte
		filePosition, ok := f.contigs[name]
		if !ok {
			// TODO: empty contig
			//log.Println(f.path + " missing contig " + name)
			f.emptyContig(done, ch)
			close(ch)
			//fmt.Println(f.path, name, "shutting down")
			return
		}

		_, err = f.rd.Seek(int64(filePosition), os.SEEK_SET)
		if err != nil {
			log.Fatal(err)
		}
		f.br.Reset(f.rd)
		for err != io.EOF {
			select {
			case <-done:
				//fmt.Println(f.path, name, "shutting down")
				close(ch)
				return
			default:
				b, err = f.br.ReadByte()
				if err == io.EOF {
					fmt.Println("FOOBAR: EOF")
					fmt.Println(f.path, name, "EOF shutting down")
					log.Fatal("FOOBAR: EOF")
					//break
				} else if err != nil {
					log.Fatal(err)
				}
				if unicode.IsSpace(rune(b)) {
					continue
				}
				if b == '>' {
					//fmt.Println(f.path, name, "New Contig")
					f.emptyContig(done, ch)
					close(ch)
					return
				}
				ch <- &Position{
					Call:       unicode.ToUpper(rune(b)),
					Coverage:   -1,
					Proportion: -1.0,
				}
			}
		}
		// TODO: Remove?
		return
	}(ch)
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

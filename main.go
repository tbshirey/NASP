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
	"unicode"
)

var refPath = flag.String("reference", "", "Path to the reference fasta against which samples are compared")
var dupPath = flag.String("duplicates", "", "Path to the duplicates file marking possible duplicated positions")

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

	// TODO: Sort analyses by Identifier
	analyses := NewSampleAnalyses(flag.Args()...)

	for _, contig := range analyses[0].(*Fasta).Contigs() {
		fmt.Println(contig)
	}

	/*
		fmt.Println(len(analyses))
		done := make(chan struct{})
		for p := range analyses[0].Contig(done, "centroid_100028") {
			fmt.Println(p)
	*/

}

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
		go func(ch chan SampleAnalysis) {
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
		}(ch)
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
}

func NewFasta(path string) *Fasta {
	fasta := &Fasta{
		path:    path,
		contigs: make(map[string]uint),
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

func (f Fasta) Contig(done chan struct{}, name string) chan *Position {
	ch := make(chan *Position)
	go func(ch chan *Position) {
		filePosition, ok := f.contigs[name]
		if !ok {
			// TODO: empty contig
			log.Fatal(f.path + " does not contain " + name)
		}
		file, err := os.Open(f.path)
		if err != nil {
			log.Fatal(err)
		}
		_, err = file.Seek(int64(filePosition), os.SEEK_SET)
		if err != nil {
			log.Fatal(err)
		}
		br := bufio.NewReader(file)
		for err != io.EOF {
			b, err := br.ReadByte()
			if err == io.EOF {
				fmt.Println("FOOBAR: EOF")
				break
			} else if err != nil {
				log.Fatal(err)
			}
			if unicode.IsSpace(rune(b)) {
				continue
			}
			if b == '>' {
				break
			}
			ch <- &Position{
				Call:       unicode.ToUpper(rune(b)),
				Coverage:   -1,
				Proportion: -1.0,
			}
		}
		// TODO: Yield empty positions
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

	file, err := os.Open(f.path)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	br := bufio.NewReader(file)
	for err != io.EOF {
		line, err = br.ReadSlice('>')
		switch err {
		default:
			log.Fatal(err)
		case bufio.ErrBufferFull, io.EOF:
			break
		case nil:
			position += uint(len(line))
			line, err = br.ReadBytes('\n')
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

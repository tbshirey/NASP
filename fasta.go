package main

import (
	"bufio"
	"bytes"
	"io"
	"log"
	"os"
	"unicode"
)

type Fasta struct {
	path string
	// TODO: aligner string
	contigs map[string]uint
	rd      *os.File
	br      *bufio.Reader
}

func NewFasta(path string, indexContigs bool) *Fasta {
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

	if indexContigs {
		fasta.index()
	}
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

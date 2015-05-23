package main

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"unicode"
)

type Reference struct {
	ref *Fasta
	dup *Fasta
}

func NewReference(refPath, dupPath string) (*Reference, error) {
	var ref, dup *Fasta
	var err error

	ref, err = NewFasta(refPath, false)
	if err != nil {
		return nil, fmt.Errorf("%s: %s", err.Error(), dupPath)
	}

	if dupPath != "" {
		// duplicates.txt is optional
		dup, err = NewFasta(dupPath, false)
		if err != nil {
			return nil, fmt.Errorf("%s: %s", err.Error(), dupPath)

		}
	}

	return &Reference{
		ref: ref,
		dup: dup,
	}, nil
}

// NextContig moves the Reader to the first position of the next contig and returns
// its name.
func (r Reference) NextContig() (name string, err error) {
	if r.dup == nil {
		return r.ref.NextContig()
	}

	// The duplicates file is optional. If present, it must be read in lockstep
	// with the reference.
	dname, derr := r.dup.NextContig()
	rname, rerr := r.ref.NextContig()

	if derr != nil || derr != rerr {
		// FIXME: Return an error warning the user the errors don't match
		// It's fine if they both return the same error such as io.EOF, but a
		// mismatch could indicate a serious error that might otherwise not be
		// reported.
		return "", err
	}

	if rname != dname {
		return "", errors.New("The duplicates file should have a corresponding contig for every contig in the reference.\n" +
			"Reference: " + rname + "\n" +
			"Duplicates: " + dname + "\n")
	}

	return rname, err
}

type Fasta struct {
	//path string
	// TODO: aligner string
	// index maps contigName -> filePosition of the first position in each contig.
	index map[string]int
	rd    *os.File
	br    *bufio.Reader
}

func NewFasta(path string, indexContigs bool) (*Fasta, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	fasta := &Fasta{
		//path:    path,
		rd: file,
		br: bufio.NewReader(file),
	}

	if indexContigs {
		// TODO: measure memory/runtime of map vs sorted list.
		fasta.index = make(map[string]int)
		fasta.indexContigs()
	}
	return fasta, nil
}

func (f Fasta) Contig(done chan struct{}, name string) chan byte {
	ch := make(chan byte, 100)
	go func(done chan struct{}, ch chan byte) {
		defer close(ch)
		var err error
		var b byte
		var r rune

		filePosition, ok := f.index[name]
		if !ok {
			return
		}

		if _, err = f.rd.Seek(int64(filePosition), os.SEEK_SET); err != nil {
			log.Fatal(err)
			//return err
		}

		f.br.Reset(f.rd)
		for err != io.EOF {
			b, err = f.br.ReadByte()
			if err != nil && err != io.EOF {
				log.Fatal(err)
				//return err
			}
			r = rune(b)
			switch {
			case unicode.IsSpace(r):
				continue
			case r == '>':
				f.br.UnreadByte()
				//return nil
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

func (f Fasta) NextContig() (name string, err error) {
	var line []byte
	for {
		line, err = f.br.ReadSlice('>')
		switch err {
		default:
			return "", err
		case bufio.ErrBufferFull:
			break
		case nil:
			if line, err = f.br.ReadBytes('\n'); err != nil {
				return "", err
			}

			//line = bytes.TrimPrefix(line, []byte("franken::"))

			if idx := bytes.IndexAny(line, " \n"); idx < 0 {
				return string(line), nil
			} else {
				return string(line[:idx]), nil
			}
		}
	}
}

/**
 * indexContigs maps the starting file position of each contig in the file.
 */
func (f Fasta) indexContigs() error {
	// TODO: Include length
	var line []byte
	var err error
	var position int

	for err != io.EOF {
		line, err = f.br.ReadSlice('>')
		switch err {
		default:
			return err
		case bufio.ErrBufferFull, io.EOF:
			break
		case nil:
			position += len(line)
			line, err = f.br.ReadBytes('\n')
			if err != nil {
				log.Fatal(err)
			}
			position += len(line)

			line = bytes.TrimPrefix(line, []byte("franken::"))

			var name string
			if idx := bytes.IndexAny(line, " \n"); idx < 0 {
				name = string(line)
			} else {
				name = string(line[:idx])
			}
			f.index[name] = position
		}
	}

	return nil
}

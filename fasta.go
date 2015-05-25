package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
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

// NextContig moves the Reader to the first position of the next contig and
// returns its name.
func (r Reference) NextContig() (string, error) {
	var rname, dname string
	var rerr, derr error

	rname, rerr = r.ref.NextContig()

	if r.dup != nil {
		// The duplicates file is optional. If present, it must be read in lockstep
		// with the reference.
		dname, derr = r.dup.NextContig()

		if derr != nil || derr != rerr {
			return "", fmt.Errorf("duplicates file: %s reference file: %s\n", derr.Error(), rerr.Error())
		}

		if rname != dname {
			return "", fmt.Errorf("The duplicates file should have a corresponding contig for every contig in the reference.\nThe following contigs were found in corresponding positions of the reference and duplicates files: `%s` and `%s`\n", rname, dname)
		}
	}

	return rname, rerr
}

func (r Reference) ReadPositions() (ref, dup []byte, isPrefix bool, rerr error) {
	var derr error

	ref, rerr = r.ref.ReadPositions()

	if r.dup != nil {
		dup, derr = r.dup.ReadPositions()

		if derr != nil || derr != rerr {
			return nil, nil, false, fmt.Errorf("duplicates file: %s reference file: %s\n", derr.Error(), rerr.Error())
		}
	}

	if rerr == bufio.ErrBufferFull {
		return ref, dup, true, nil
	}
	return ref, dup, false, rerr
}

type Fasta struct {
	//path string
	// TODO: aligner string
	// index maps contigName -> filePosition of the first position in each contig.
	name  string
	index map[string]int64
	rd    *os.File
	br    *bufio.Reader
	buf   *bytes.Buffer
}

func NewFasta(path string, indexContigs bool) (*Fasta, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	fasta := &Fasta{
		//path:    path,
		name: filepath.Base(path),
		rd:   file,
		br:   bufio.NewReader(file),
		buf:  &bytes.Buffer{},
	}

	if indexContigs {
		// TODO: measure memory/runtime of map vs sorted list.
		fasta.index = make(map[string]int64)
		fasta.indexContigs()
	}
	return fasta, nil
}

/*
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
*/

func (f Fasta) Name() string {
	return f.name
}

/**
 * NextContig advances the io.Reader to the first position of the next contig
 * and returns its name.
 *
 * TODO: document name excludes description and prefix
 * >frankenfasta::name description
 * gatc...
 */
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

func (f Fasta) ReadPositions() ([]byte, error) {
	return f.br.ReadSlice('>')
}

func (f Fasta) SeekContig(name string) error {
	if filePosition, ok := f.index[name]; ok {
		_, err := f.rd.Seek(filePosition, os.SEEK_SET)
		return err
	}

	_, err := f.rd.Seek(0, os.SEEK_END)
	return err
}

/*
func (f Fasta) ScanPositions(name string) ([]byte, error) {
	filePosition, ok := f.index[name]
	if !ok {
		return nil, nil
	}
	fmt.Printf("Scanning %s at %d\n", name, filePosition)
	f.rd.Seek(filePosition, os.SEEK_SET)
	f.br.Reset(f.rd)
	return f.br.ReadBytes('>')
}
*/

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
			position += len(line)
			break
		case nil:
			position += len(line)
			// NOTE: This saves a few megabytes of RAM over ReadBytes,
			// but could fail with a ErrBufferFull
			line, err = f.br.ReadSlice('\n')
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
			f.index[name] = int64(position)
		}
	}

	return nil
}

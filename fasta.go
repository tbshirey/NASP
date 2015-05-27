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

func (r Reference) ReadPositions(n int) (ref, dup []byte, isPrefix bool, rerr error) {
	var derr error

	ref, rerr = r.ref.ReadPositions(n)

	if r.dup != nil {
		dup, derr = r.dup.ReadPositions(n)

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
	name     string
	index    map[string]int64
	rd       *os.File
	br       *bufio.Reader
	isPrefix bool
	buf      *bytes.Buffer
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
func (f *Fasta) NextContig() (name string, err error) {
	f.isPrefix = true
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

// ReadPositions returns the next `n` positions or nil when the contig is
// exhausted. `n` can be any arbitrary buffer size.
func (f Fasta) ReadPositions(n int) ([]byte, error) {
	if f.isPrefix {
		for f.buf.Len() < n {
			line, _, err := f.br.ReadLine()
			if err != nil {
				line = nil
				break
			}
			if peek, _ := f.br.Peek(1); peek[0] == '>' {
				f.isPrefix = false
				break
			}
			// err is always nil
			f.buf.Write(bytes.ToUpper(line))
		}
	}
	return f.buf.Next(n), nil
}

// SeekContig moves the
func (f *Fasta) SeekContig(name string) (err error) {
	if filePosition, ok := f.index[name]; ok {
		// ReadPositions will yield values from the new contig.
		f.isPrefix = true
		_, err = f.rd.Seek(filePosition, os.SEEK_SET)
		// Purge any old data from the previous file position
		// TODO: The contigs will probably be typically sequential, in which
		// case it might be efficient to advance the buffers instead of purging
		// them.
		f.br.Reset(f.rd)
		f.buf.Reset()
	} else {
		f.isPrefix = false
	}
	return err
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

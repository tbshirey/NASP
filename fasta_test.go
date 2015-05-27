package main

import (
	"bytes"
	"io"
	"io/ioutil"
	"os"
	"testing"
)

func tmpFasta(t *testing.T) *os.File {
	file, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := file.Write([]byte(">ContigA\nGATC\nABCD\r\nefgh\n>ContigB\nCTAGDCBA")); err != nil {
		t.Fatal(err)
	}
	return file
}

func closeFasta(file *os.File) {
	file.Close()
	os.Remove(file.Name())
}

func TestNewReferenceInvalidReference(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	_, err := NewReference("invalid", "")
	if err == nil {
		t.Fatal(err)
	}

	t.Log(err.Error())
}

func TestNewReferenceInvalidDuplicates(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	dupFasta := tmpFasta(t)
	defer closeFasta(dupFasta)

	_, err := NewReference(refFasta.Name(), "invalid")
	if err == nil {
		t.Fatal(err)
	}

	t.Log(err.Error())
}

func TestNewReferenceOnly(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	if ref, err := NewReference(refFasta.Name(), ""); ref == nil || err != nil {
		t.Fatal(err)
	}
}

func TestNewReferenceWithDuplicates(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	dupFasta := tmpFasta(t)
	defer closeFasta(dupFasta)

	if ref, err := NewReference(refFasta.Name(), dupFasta.Name()); ref == nil || err != nil {
		t.Fatal(err)
	}
}

func TestNext(t *testing.T) {
	file := tmpFasta(t)
	defer closeFasta(file)

	fasta, err := NewFasta(file.Name(), false)
	if err != nil {
		t.Fatal(err)
	}

	expect := "ContigA"
	name, err := fasta.NextContig()
	if name != expect || err != nil {
		t.Fatalf("Expected: %s, nil Observed: %s, %s", name, err)
	}

	expect = "ContigB"
	name, err = fasta.NextContig()
	if name != expect || err != nil {
		t.Fatalf("Expected: %s, nil Observed: %s, %s", name, err)
	}

	expect = ""
	name, err = fasta.NextContig()
	if name != expect || err != io.EOF {
		t.Fatalf("Expected: %s, io.EOF Observed: %s, %s", name, err)
	}
}

func TestIndex(t *testing.T) {
	file, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		file.Close()
		os.Remove(file.Name())
	}()

	if _, err := file.Write([]byte(">ContigA\nGATC\nABCD\r\nefgh\r\n>ContigB\nCTAGDCBA")); err != nil {
		t.Fatal(err)
	}

	fasta, err := NewFasta(file.Name(), true)
	if err != nil {
		t.Fatal(err)
	}

	done := make(chan struct{})
	defer close(done)

	contigName := "ContigA"
	if err = fasta.SeekContig(contigName); err != nil {
		t.Fatal(err)
	}
	expect := []byte("GATCABCDEFGH")
	actual, err := fasta.ReadPositions(4)
	if !bytes.Equal(expect, actual) {
		t.Fatalf("%s expected: %s actual: %s", contigName, expect, actual)
	}
	// TODO: Expected it yields empty positions when the contig is exhausted

	contigName = "ContigB"
	if err = fasta.SeekContig(contigName); err != nil {
		t.Fatal(err)
	}
	expect = []byte("CTAGDCBA")
	actual, err = fasta.ReadPositions(4)
	if !bytes.Equal(expect, actual) {
		t.Fatalf("%s expected: %c actual: %c", contigName, expect, actual)
	}

	// TODO: Expected it yields empty positions when the file is exhausted
	// until the channel is closed
}

package main

import (
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"testing"
	"time"
)

func TestNext(t *testing.T) {
	file, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		file.Close()
		os.Remove(file.Name())
	}()

	if _, err := file.Write([]byte(">ContigA\nGATC\nABCD\n>ContigB\nCTAGDCBA")); err != nil {
		t.Fatal(err)
	}

	fasta, err := NewFasta(file.Name(), false)
	if err != nil {
		t.Fatal(err)
	}

	expect := "ContigA"
	name, err := fasta.Next()
	if name != expect || err != nil {
		t.Fatalf("Expected: %s, nil Observed: %s, %s", name, err)
	}

	expect = "ContigB"
	name, err = fasta.Next()
	if name != expect || err != nil {
		t.Fatalf("Expected: %s, nil Observed: %s, %s", name, err)
	}

	expect = ""
	name, err = fasta.Next()
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

	if _, err := file.Write([]byte(">ContigA\nGATC\nABCD\n>ContigB\nCTAGDCBA")); err != nil {
		t.Fatal(err)
	}

	fasta, err := NewFasta(file.Name(), true)
	if err != nil {
		t.Fatal(err)
	}

	done := make(chan struct{})
	defer close(done)

	contigName := "ContigA"
	c := fasta.Contig(done, contigName)
	for i, expect := range "GATCABCD" {
		select {
		case <-time.After(time.Second):
			t.Fatalf("timeout: Expected %s position %d expected: %c", contigName, i+1, expect)
		case actual := <-c:
			if expect != rune(actual) {
				t.Fatalf("Expected %s position %d expected: %c actual: %c", contigName, i+1, expect, actual)
			}
		}
	}

	// TODO: Expected it yields empty positions when the contig is exhausted

	contigName = "ContigB"
	c = fasta.Contig(done, contigName)
	for i, expect := range "CTAGDCBA" {
		select {
		case <-time.After(time.Second):
			t.Fatalf("timeout: Expected %s position %d expected: %c", contigName, i+1, expect)
		case actual := <-c:
			if expect != rune(actual) {
				t.Fatalf("Expected %s position %d expected: %c actual: %c", contigName, i+1, expect, actual)
			}
		}
	}

	// TODO: Expected it yields empty positions when the file is exhausted
	// until the channel is closed
	fmt.Println(<-c)
}

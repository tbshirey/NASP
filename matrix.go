package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"os/signal"
	"strconv"
	"sync"
)

func master(wg sync.WaitGroup, ch chan []byte) {
	defer wg.Done()

	file, err := os.Create("master.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	buf := bufio.NewWriter(file)
	defer buf.Flush()

	for row := range ch {
		buf.Write(row)
	}
}

type LineBuilder struct {
	master      bytes.Buffer
	missingdata bytes.Buffer
}

func (l LineBuilder) Dispatch(chunk chan Chunk) {

}

func lineDispatcher(chunks chan Chunk) {
	var prevContigName string
	var positionNum int64

	for chunk := range chunks {
		positions := <-chunk.positionsChan
		for _, position := range positions {
			if prevContigName != chunk.contigName {
				positionNum = 0
			}
			positionNum++

			select {
			default:
				//case masterCh <- buildMasterRow():
				//case allRefCh <- buildAllRefCh():
			}
			positionPool.Put(position)
		}
		positionsPool.Put(positions)
	}
}

func masterMatrixLine(position *Position) {

}

// Put 1 positionsPool, maxBufSize positionPool
func writeMaster(chunks chan Chunk, numSamples int) {
	interruptChan := make(chan os.Signal, 1)
	signal.Notify(interruptChan, os.Interrupt)

	file, err := os.Create("master.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	//buf := bytes.Buffer{}
	buf := bufio.NewWriter(file)
	defer buf.Flush()

	if _, err := buf.Write([]byte("LocusID\tReference\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\n")); err != nil {
		log.Fatal(err)
	}

	numberColumns := make([]byte, 0, 100)

	var prevContigName string
	var positionNum int64

	for {
		select {
		case chunk, ok := <-chunks:
			if !ok {
				return
			}

			fmt.Println("Pending", len(chunks))
			positions := <-chunk.positionsChan

			for _, position := range positions {

				// Reset the position counter when starting a new contig.
				if chunk.contigName != prevContigName {
					positionNum = 0
					prevContigName = chunk.contigName
				}
				positionNum++

				// LocusID
				buf.WriteString(chunk.contigName)
				numberColumns = numberColumns[:0]
				numberColumns = append(numberColumns, ':', ':')
				numberColumns = strconv.AppendInt(numberColumns, positionNum, 10)
				buf.Write(numberColumns)
				buf.WriteByte('\t')

				for _, call := range position.callStr {
					buf.WriteByte(call)
					buf.WriteByte('\t')
				}

				numberColumns = numberColumns[:0]
				numberColumns = strconv.AppendInt(numberColumns, int64(position.calledSnp), 10)
				numberColumns = append(numberColumns, '\t')
				// #Indelcall
				numberColumns = strconv.AppendInt(numberColumns, int64(position.calledReference), 10)
				numberColumns = append(numberColumns, '\t')
				numberColumns = strconv.AppendInt(numberColumns, int64(position.a), 10)
				numberColumns = append(numberColumns, '\t')
				numberColumns = strconv.AppendInt(numberColumns, int64(position.c), 10)
				numberColumns = append(numberColumns, '\t')
				numberColumns = strconv.AppendInt(numberColumns, int64(position.g), 10)
				numberColumns = append(numberColumns, '\t')
				numberColumns = strconv.AppendInt(numberColumns, int64(position.t), 10)
				numberColumns = append(numberColumns, '\t')
				numberColumns = strconv.AppendInt(numberColumns, int64(position.n), 10)
				numberColumns = append(numberColumns, '\t')
				buf.Write(numberColumns)

				buf.WriteString(chunk.contigName)
				buf.WriteByte('\t')

				numberColumns = numberColumns[:0]
				numberColumns = strconv.AppendInt(numberColumns, int64(positionNum), 10)
				numberColumns = append(numberColumns, '\t')

				buf.Write(position.callWasMade)
				buf.WriteByte('\t')
				buf.Write(position.passedDepthFilter)
				buf.WriteByte('\t')
				buf.Write(position.passedProportionFilter)
				buf.WriteByte('\n')

				positionPool.Put(position)
			}
			positionsPool.Put(positions)
		case <-interruptChan:
			return
		}
	}
}

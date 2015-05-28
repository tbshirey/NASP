package main

import (
	"bufio"
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

func lineDispatcher(ch chan chan []*Position) {
	for positions := range ch {
		for _, position := range <-positions {
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

				if chunk.contigName != prevContigName {
					positionNum = 0
				}

				positionNum++

				buf.WriteString(chunk.contigName)

				numberColumns = numberColumns[:0]
				numberColumns = append(numberColumns, ':', ':')
				strconv.AppendInt(numberColumns, positionNum, 10)
				buf.Write(numberColumns)
				buf.WriteByte('\t')

				for _, call := range position.callStr {
					buf.WriteByte(call)
					buf.WriteByte('\t')
				}

				numberColumns = numberColumns[:0]
				strconv.AppendInt(numberColumns, int64(position.calledSnp), 10)
				numberColumns = append(numberColumns, '\t')
				// #Indelcall
				strconv.AppendInt(numberColumns, int64(position.calledReference), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.a), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.c), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.g), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.t), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.n), 10)
				numberColumns = append(numberColumns, '\t')
				buf.Write(numberColumns)

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

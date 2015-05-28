package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"os/signal"
	"strconv"
)

var (
	NEWLINE = []byte("\n")
	TAB     = []byte("\t")
)

func writeMaster(ch chan chan []*Position, numSamples int) {
	/*
				'LocusID': "{0}::{1}".format(contig_name, position),
		        #         'Reference': row.call_str[0],
		        #         '#SNPcall': row.called_snp,
		        #         # TODO: replace with n/a
		        #         '#Indelcall': '0',
		        #         '#Refcall': row.called_reference,
		        #         '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
		        #         '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
		        #         '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
		        #         '#A': row.num_A,
		        #         '#C': row.num_C,
		        #         '#G': row.num_G,
		        #         '#T': row.num_T,
		        #         # TODO: replace with n/a
		        #         '#Indel': '0',
		        #         '#NXdegen': row.num_N,
		        #         'Contig': contig_name,
		        #         'Position': position,
		        #         'InDupRegion': row.is_reference_duplicated,
		        #         'SampleConsensus': row.is_all_passed_consensus,
		        #         'CallWasMade': row.CallWasMade,
		        #         'PassedDepthFilter': row.PassedDepthFilter,
		        #         'PassedProportionFilter': row.PassedProportionFilter,
		        #         'Pattern': "".join(row.Pattern)
	*/
	interruptChan := make(chan os.Signal, 1)
	signal.Notify(interruptChan, os.Interrupt)

	file, err := os.Create("master.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	// Compressed Writer
	/*
		cw := gzip.NewWriter(file)
		defer func() {
			log.Println("master.tsv.gz goroutine shutting down")
			if err := cw.Close(); err != nil {
				log.Fatal(err)
			}
			if err = file.Close(); err != nil {
				log.Fatal(err)
			}
		}()
	*/

	//buf := bytes.Buffer{}
	buf := bufio.NewWriter(file)
	defer buf.Flush()

	if _, err := buf.Write([]byte("LocusID\tReference\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\n")); err != nil {
		log.Fatal(err)
	}
	numberColumns := make([]byte, 0, 100)
	for {
		select {
		case positions, ok := <-ch:
			if !ok {
				return
			}
			fmt.Println("Pending", len(ch))
			p := <-positions
			for _, position := range p {
				/*
					buf.Reset()
					buf.Write(row.Contig)
					if _, err = buf.WriteTo(cw); err != nil {
						log.Fatal(err)
					}
				*/

				for _, call := range position.callStr {
					buf.WriteByte(call)
					buf.WriteByte('\t')
				}

				numberColumns = numberColumns[:0]
				strconv.AppendInt(numberColumns, int64(position.a), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.c), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.g), 10)
				numberColumns = append(numberColumns, '\t')
				strconv.AppendInt(numberColumns, int64(position.t), 10)
				numberColumns = append(numberColumns, '\t')

				buf.Write(numberColumns)

				//buf.Write(position.callStr)
				//fmt.Fprintf(buf, "%d\t%d\t%d\t%d", position.a, position.c, position.g, position.g)
				buf.Write(position.callWasMade)
				//buf.Write(TAB)
				buf.WriteByte('\t')
				buf.Write(position.passedDepthFilter)
				buf.WriteByte('\t')
				//	buf.Write(TAB)
				buf.Write(position.passedProportionFilter)
				buf.WriteByte('\n')
				//buf.Write(NEWLINE)
				//buf.WriteTo(cw)
				//buf.Reset()

				/*
					if _, err = cw.Write(position.callStr); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(position.callWasMade); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(TAB); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(position.passedDepthFilter); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(TAB); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(position.passedProportionFilter); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(TAB); err != nil {
						log.Fatal(err)
					}
					if _, err = cw.Write(NEWLINE); err != nil {
						log.Fatal(err)
					}
				*/

				//buf.Reset()
				positionPool.Put(position)
			}
			positionsPool.Put(p)
		case <-interruptChan:
			return
		}
	}
}

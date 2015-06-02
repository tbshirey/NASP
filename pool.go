package main

import (
	"bytes"
	"sync"
)

var callsPool = sync.Pool{
	New: func() interface{} {
		calls := make([][]byte, NUM_SAMPLES)
		for i := range calls {
			calls[i] = make([]byte, defaultBufSize)
		}
		return calls
	},
}

var positionsPool = sync.Pool{
	New: func() interface{} {
		return make([]*Position, defaultBufSize)
	},
}

type Position struct {
	// General Stats
	isAllCalled           bool
	isReferenceClean      bool
	isReferenceDuplicated bool
	isAllPassCoverage     bool
	isAllPassProportion   bool
	isAllPassConsensus    bool
	isAllQualityBreadth   bool
	isBestSnp             bool
	isMissingMatrix       bool

	//all_sample_stats=all_sample_stats,

	// Missing Data Matrix condition - at least one SampleAnalysis passes quality_breadth and is a SNP.
	//is_missing_matrix=is_missing_matrix,

	// NASP Master Matrix
	// Counters
	wasCalled        int64
	calledReference  int64
	calledSnp        int64
	calledDegen      int64
	passedCoverage   int64
	passedProportion int64
	a                int64
	c                int64
	g                int64
	t                int64
	n                int64

	// Strings
	callStr []byte
	// TODO
	//maskedCallStr          []byte
	callWasMade            []byte
	passedDepthFilter      []byte
	passedProportionFilter []byte
	pattern                []byte
}

var positionPool = sync.Pool{
	New: func() interface{} {
		return &Position{
			isAllCalled: true,
			//isReferenceClean:      true,
			//isReferenceDuplicated: true,
			isAllPassCoverage:   true,
			isAllPassProportion: true,
			isAllPassConsensus:  true,
			isAllQualityBreadth: true,
			isBestSnp:           true,
			//isMissingMatrix: true,

			callStr: make([]byte, NUM_SAMPLES+1),
			// TODO
			//maskedCallStr:          make([]byte, numSamples),
			callWasMade:            make([]byte, NUM_SAMPLES),
			passedDepthFilter:      make([]byte, NUM_SAMPLES),
			passedProportionFilter: make([]byte, NUM_SAMPLES),
			pattern:                make([]byte, NUM_SAMPLES+1),
		}
	},
}

var statsPool = sync.Pool{
	New: func() interface{} {
		return make(SampleStats, NUM_SAMPLES)
	},
}

var linePool = sync.Pool{
	New: func() interface{} {
		return &bytes.Buffer{}
	},
}

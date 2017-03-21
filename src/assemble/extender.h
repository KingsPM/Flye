//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <deque>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "chimera.h"

struct ContigPath
{
	ContigPath(): 
		circular(false) {}

	std::vector<FastaRecord::Id> reads;
	bool circular;
};

class Extender
{
public:
	Extender(const SequenceContainer& readsContainer, 
			 OverlapContainer& ovlpContainer,
			 int coverage, int genomeSize):
		_readsContainer(readsContainer), 
		_ovlpContainer(ovlpContainer),
		_chimDetector(coverage, readsContainer, ovlpContainer),
		_coverage(coverage), _genomeSize(genomeSize),
		_progress(genomeSize)
	{}

	void assembleContigs();
	const std::vector<ContigPath>& getContigPaths() const
		{return _contigPaths;}

private:
	const SequenceContainer& _readsContainer;
	OverlapContainer& _ovlpContainer;
	ChimeraDetector   _chimDetector;
	const int 		  _coverage;
	const int 		  _genomeSize;
	ProgressPercent   _progress;

	FastaRecord::Id stepRight(FastaRecord::Id readId);
	ContigPath extendContig(FastaRecord::Id startingRead);

	int   countRightExtensions(FastaRecord::Id readId);
	bool  extendsRight(const OverlapRange& ovlp);


	std::vector<ContigPath> 				 _contigPaths;
	std::unordered_set<FastaRecord::Id>		 _visitedReads;
	bool _rightExtension;
	int  _assembledSequence;
};

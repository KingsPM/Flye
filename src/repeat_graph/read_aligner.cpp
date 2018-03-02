//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "read_aligner.h"
#include "../common/parallel.h"

namespace
{
	struct Chain
	{
		Chain(): score(0) {}
		std::vector<const EdgeAlignment*> aln;
		int32_t score;
	};
}

//Give alignments to separate edges for a single read, merges them
//into non-overlapping chains (could be more than one chain per read
//in case of chimera) with maximum score
std::vector<GraphAlignment>
	ReadAligner::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 	 const std::vector<EdgeAlignment>& ovlps) const
{
	static const int32_t MAX_DISCORDANCE = 
		std::max(Config::get("maximum_jump") / Config::get("far_jump_rate"),
				 Config::get("max_separation"));
	static const int32_t MAX_JUMP = Config::get("maximum_jump");
	static const int32_t ALN_GAP = Config::get("read_align_gap");
	static const int32_t PENALTY_WND = Config::get("penalty_window");

	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxScore = 0;
		Chain* maxChain = nullptr;
		for (auto& chain : activeChains)
		{
			const OverlapRange& nextOvlp = edgeAlignment.overlap;
			const OverlapRange& prevOvlp = chain.aln.back()->overlap;

			int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
			int32_t graphDiff = nextOvlp.extBegin +
								prevOvlp.extLen - prevOvlp.extEnd;

			if (chain.aln.back()->edge->nodeRight == edgeAlignment.edge->nodeLeft &&
				MAX_JUMP > readDiff && readDiff > 0 &&
				MAX_JUMP > graphDiff && graphDiff > 0  &&
				abs(readDiff - graphDiff) < MAX_DISCORDANCE)
			{
				int32_t gapScore = -(readDiff - ALN_GAP) / PENALTY_WND;
				if (readDiff < ALN_GAP) gapScore = 1;
				//if (chain.aln.back()->segment.end != 
				//	edgeAlignment.segment.start) gapScore -= 10;
				//int32_t ovlpScore = !edgeAlignment.edge->isLooped() ? nextOvlp.score : 10;
				int32_t score = chain.score + nextOvlp.score + gapScore;
				if (score > maxScore)
				{
					maxScore = score;
					maxChain = &chain;
				}
			}
		}
		
		if (maxChain)
		{
			newChains.push_back(*maxChain);
			maxChain->aln.push_back(&edgeAlignment);
			maxChain->score = maxScore;
		}

		activeChains.splice(activeChains.end(), newChains);
		activeChains.push_back(Chain());
		activeChains.back().aln.push_back(&edgeAlignment);
		activeChains.back().score = edgeAlignment.overlap.score;
	}

	//choosing optimal(ish) set of alignments
	std::vector<GraphAlignment> acceptedAlignments;
	std::vector<Chain> sortedChains(activeChains.begin(), activeChains.end());
	std::sort(sortedChains.begin(), sortedChains.end(),
			  [](const Chain& c1, const Chain& c2)
			  {return c1.score > c2.score;});
	for (auto& chain : sortedChains)
	{
		int32_t alnLen = chain.aln.back()->overlap.curEnd - 
					 	 chain.aln.front()->overlap.curBegin;
		if (alnLen < Parameters::get().minimumOverlap) continue;

		//check if it overlaps with other accepted chains
		bool overlaps = false;
		for (auto& existAln : acceptedAlignments)
		{
			int32_t existStart = existAln.front().overlap.curBegin;
			int32_t existEnd = existAln.back().overlap.curEnd;
			int32_t curStart = chain.aln.front()->overlap.curBegin;
			int32_t curEnd = chain.aln.back()->overlap.curEnd;

			int32_t overlapRate = std::min(curEnd, existEnd) - 
									std::max(curStart, existStart);
			if (overlapRate > (int)Config::get("max_separation")) overlaps = true;
		}
		if (!overlaps) 
		{
			acceptedAlignments.emplace_back();
			for (auto& aln : chain.aln) acceptedAlignments.back().push_back(*aln);
		}
	}

	return acceptedAlignments;
}

void ReadAligner::alignReads()
{
	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, SequenceSegment>> idToSegment;
	SequenceContainer pathsContainer;

	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand()) continue;

		//add each edge sequence variant to the database
		for (auto& segment : edge->seqSegments)
		{
			size_t len = segment.end - segment.start;
			auto sequence = _asmSeqs.getSeq(segment.seqId)
										.substr(segment.start, len);
			auto& newRec = pathsContainer.addSequence(sequence, "");

			idToSegment[newRec.id] = {edge, segment};
			idToSegment[newRec.id.rc()] = {_graph.complementEdge(edge), 
										   segment.complement()};
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(pathsContainer, 
						   (int)Config::get("read_align_kmer_sample"));
	pathsIndex.countKmers(1, /* genome size*/ 0);
	pathsIndex.buildIndex(1, (int)Config::get("read_align_max_kmer"));
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									(int)Config::get("maximum_jump"),
									(int)Config::get("max_separation"),
									/*no overhang*/0,
									(int)Config::get("read_align_gap"),
						 			(float)Config::get("read_align_max_divergence"),
									/*keep alignment*/ false);

	//readsOverlapper.setDumpFile("../graph-reads-div.txt");
	OverlapContainer readsOverlaps(readsOverlapper, _readSeqs, 
								   /*onlyMax*/ false);

	std::vector<FastaRecord::Id> allQueries;
	int64_t totalLength = 0;
	for (auto& readId : _readSeqs.getIndex())
	{
		//if (_readSeqs.seqLen(readId.first) > Constants::maxSeparation)
		if (_readSeqs.seqLen(readId.first) > Parameters::get().minimumOverlap)
		{
			totalLength += _readSeqs.seqLen(readId.first);
			allQueries.push_back(readId.first);
		}
	}
	std::mutex indexMutex;
	int numAligned = 0;
	int64_t alignedLength = 0;
	std::function<void(const FastaRecord::Id&)> alignRead = 
	[this, &indexMutex, &numAligned, &readsOverlaps, 
		&idToSegment, &pathsContainer, &alignedLength] 
	(const FastaRecord::Id& seqId)
	{
		auto overlaps = readsOverlaps.seqOverlaps(seqId);
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : overlaps)
		{
			alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
								  idToSegment[ovlp.extId].second});
		}
		std::sort(alignments.begin(), alignments.end(),
		  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			{return e1.overlap.curBegin < e2.overlap.curBegin;});
		auto readChains = this->chainReadAlignments(pathsContainer, alignments);

		if (readChains.empty()) return;
		indexMutex.lock();
		++numAligned;
		for (auto& chain : readChains) 
		{
			_readAlignments.push_back(chain);
			alignedLength += chain.back().overlap.curEnd - 
							 chain.front().overlap.curBegin;
		}
		indexMutex.unlock();
	};

	processInParallel(allQueries, alignRead, 
					  Parameters::get().numThreads, true);

	/*for (auto& aln : _readAlignments)
	{
		if (aln.size() > 1)
		{
			std::string alnStr;
			int switches = 0;
			for (size_t i = 0; i < aln.size() - 1; ++i)
			{
				if (aln[i].segment.end != aln[i + 1].segment.start) ++switches;
			}

			for (auto& edge : aln)
			{
				alnStr += std::to_string(edge.edge->edgeId.signedId()) + " ("
					+ std::to_string(edge.overlap.curRange()) + " " 
					+ std::to_string(edge.overlap.score) + " " 
					+ std::to_string((float)edge.overlap.score / 
						edge.overlap.curRange()) + ") -> ";
			}
			alnStr.erase(alnStr.size() - 4);
			alnStr += " s: " + std::to_string(switches);
			FastaRecord::Id readId = aln.front().overlap.curId;
			Logger::get().debug() << "Aln " << _readSeqs.seqName(readId);
			Logger::get().debug() << "\t" << alnStr;
		}
	}*/

	Logger::get().debug() << "Aligned " << numAligned << " / " 
		<< allQueries.size();
	Logger::get().debug() << "Aligned length " << alignedLength << " / " 
		<< totalLength << " " << (float)alignedLength / totalLength;
}

//updates alignments with respect to the new graph
void ReadAligner::updateAlignments()
{
	std::vector<GraphAlignment> newAlignments;
	for (auto& aln : _readAlignments)
	{
		GraphAlignment curAlignment;
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			if (!_graph.hasEdge(aln[i].edge)) continue;

			curAlignment.push_back(aln[i]);
			if (!_graph.hasEdge(aln[i + 1].edge) ||
				aln[i].edge->nodeRight != aln[i + 1].edge->nodeLeft)
			{
				newAlignments.push_back(curAlignment);
				curAlignment.clear();
			}
		}

		if (_graph.hasEdge(aln.back().edge)) curAlignment.push_back(aln.back());
		if (!curAlignment.empty()) newAlignments.push_back(curAlignment);
	}

	_readAlignments = newAlignments;
}

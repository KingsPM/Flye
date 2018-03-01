//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"

#include "../ksw2/ksw2.h"

namespace
{
	float kswAlign(const std::vector<uint8_t>& tseq, 
				   const std::vector<uint8_t>& qseq, int matchScore, 
			   	   int misScore, int gapOpen, int gapExtend, int bandWidth)
	{
		//substitution matrix
		int8_t a = matchScore;
		int8_t b = misScore < 0 ? misScore : -misScore; // a > 0 and b < 0
		int8_t subsMat[] = {a, b, b, b, 0, 
						  	b, a, b, b, 0, 
						  	b, b, a, b, 0, 
						  	b, b, b, a, 0, 
						  	0, 0, 0, 0, 0};

		//encode as nucleotide ids
		/*uint8_t nucToId[256];
		memset(nucToId, 4, 256);
		nucToId[(uint8_t)'A'] = nucToId[(uint8_t)'a'] = 0; 
		nucToId[(uint8_t)'C'] = nucToId[(uint8_t)'c'] = 1;
		nucToId[(uint8_t)'G'] = nucToId[(uint8_t)'g'] = 2; 
		nucToId[(uint8_t)'T'] = nucToId[(uint8_t)'t'] = 3;*/

		//uint8_t *ts, *qs;
		//ts = (uint8_t*)malloc(tseq.size());
		//qs = (uint8_t*)malloc(qseq.size());
		//for (size_t i = 0; i < tseq.size(); ++i) ts[i] = nucToId[(uint8_t)tseq[i]];
		//for (size_t i = 0; i < qseq.size(); ++i) qs[i] = nucToId[(uint8_t)qseq[i]];

		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ksw_extz_t));
		const int NUM_NUCL = 5;
		const int Z_DROP = -1;
		const int FLAG = KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
		//ksw_extf2_sse(0, qseq.size(), &qseq[0], tseq.size(), &tseq[0], matchScore,
		//		 	  misScore, gapOpen, bandWidth, Z_DROP, &ez);
		ksw_extz2_sse(0, qseq.size(), &qseq[0], tseq.size(), &tseq[0], NUM_NUCL,
				 	  subsMat, gapOpen, gapExtend, bandWidth, Z_DROP, FLAG, &ez);
		
		//decode CIGAR
		/*std::string strQ;
		for (auto x : qseq) strQ += "ACGT"[x];
		std::string strT;
		for (auto x : tseq) strT += "ACGT"[x];
		std::string alnQry;
		std::string alnTrg;*/

		size_t posQry = 0;
		size_t posTrg = 0;
        int matches = 0;
		int alnLength = 0;
		for (size_t i = 0; i < (size_t)ez.n_cigar; ++i)
		{
			int size = ez.cigar[i] >> 4;
			char op = "MID"[ez.cigar[i] & 0xf];
			alnLength += size;

        	if (op == 'M')
			{
				for (size_t i = 0; i < (size_t)size; ++i)
				{
					if (tseq[posTrg + i] == qseq[posQry + i]) ++matches;
				}

                //alnQry += strQ.substr(posQry, size);
                //alnTrg += strT.substr(posTrg, size);
				posQry += size;
				posTrg += size;
			}
            else if (op == 'I')
			{
                //alnQry += strQ.substr(posQry, size);
                //alnTrg += std::string(size, '-');
                posQry += size;
			}
            else //D
			{
                //alnQry += std::string(size, '-');
                //alnTrg += strT.substr(posTrg, size);
				posTrg += size;
			}
		}
		/*for (size_t i = 0; i < alnTrg.size() / 100 + 1; ++i)
		{
			std::cout << alnTrg.substr(i * 100, 100) << "\n" 
				<< alnQry.substr(i * 100, 100) << "\n\n";
		}*/
        float errRate = 1 - float(matches) / alnLength;

		free(ez.cigar);
		return errRate;
	}
}


//pre-filtering
bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
								int32_t curLen, int32_t extLen) const
{	
	if (_checkOverhang && 
		std::min(curPos, extPos) > _maxOverhang) return false;

	if (extPos > extLen - _minOverlap ||
	    curPos > curLen - _minOverlap) return false;

	return true;
}


OverlapDetector::JumpRes 
OverlapDetector::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext) const
{
	static const int CLOSE_JUMP = Config::get("close_jump_rate");
	static const int FAR_JUMP = Config::get("far_jump_rate");

	if (curNext - curPrev > _maxJump) return J_END;

	if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
		0 < extNext - extPrev && extNext - extPrev < _maxJump)
	{
		float jumpDivergence = abs((curNext - curPrev) - (extNext - extPrev));
		float jumpLength = ((curNext - curPrev) + (extNext - extPrev)) / 2;

		if (jumpDivergence < _maxJump / CLOSE_JUMP &&
			jumpLength < _gapSize)
		{
			return J_CLOSE;
		}
		if (jumpDivergence < _maxJump / FAR_JUMP)
		{
			return J_FAR;
		}
	}
	return J_INCONS;
}


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp) const
{
	static const int OVLP_DIVERGENCE = Config::get("overlap_divergence_rate");
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	float ovlpDiff = abs(ovlp.curRange() - ovlp.extRange());
	float ovlpLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (ovlpDiff > ovlpLength / OVLP_DIVERGENCE)
	{
		return false;
	}

	if (_checkOverhang && ovlp.curId != ovlp.extId.rc())
	{
		if (std::min(ovlp.curBegin, ovlp.extBegin) > 
			_maxOverhang) 
		{
			return false;
		}
		if (std::min(ovlp.curLen - ovlp.curEnd, ovlp.extLen - ovlp.extEnd) > 
			_maxOverhang)
		{
			return false;
		}
	}

	
	return true;
}

namespace
{
	struct DPRecord
	{
		DPRecord() {}
		DPRecord(const OverlapRange& ovlp): ovlp(ovlp) {}

		OverlapRange ovlp;
		std::vector<int32_t> shifts;
	};
}

std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool uniqueExtensions) const
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> activePaths;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> completedPaths;
	std::set<size_t> eraseMarks;
	int32_t curLen = fastaRec.sequence.length();

	//for all kmers in this read
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		int32_t curPos = curKmerPos.position;

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if (extReadPos.readId == fastaRec.id &&
				extReadPos.position == curPos) continue;

			int32_t extLen = _seqContainer.seqLen(extReadPos.readId);
			if (extLen < (int32_t)_minOverlap) continue;

			int32_t extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];

			size_t maxCloseId = 0;
			size_t maxFarId = 0;
			int32_t maxCloseScore = 0;
			int32_t maxFarScore = 0;
			bool extendsClose = false;
			bool extendsFar = false;
			eraseMarks.clear();

			//searching for longest possible extension
			for (size_t pathId = 0; pathId < extPaths.size(); ++pathId)
			{
				JumpRes jumpResult = 
					this->jumpTest(extPaths[pathId].ovlp.curEnd, curPos,
								   extPaths[pathId].ovlp.extEnd, extPos);

				static const int PENALTY_WND = Config::get("penalty_window");
				int32_t jumpLength = curPos - extPaths[pathId].ovlp.curEnd;
				int32_t gapScore = -(jumpLength - _gapSize) / PENALTY_WND;
				if (jumpLength < _gapSize) gapScore = 1;

				int32_t jumpScore = extPaths[pathId].ovlp.score + gapScore;

				switch (jumpResult)
				{
					case J_END:
						eraseMarks.insert(pathId);
						if (this->overlapTest(extPaths[pathId].ovlp))
						{
							completedPaths[extReadPos.readId]
									.push_back(extPaths[pathId]);
						}
						break;
					case J_INCONS:
						break;
					case J_CLOSE:
						eraseMarks.insert(pathId);
						if (jumpScore > maxCloseScore)
						{
							extendsClose = true;
							maxCloseId = pathId;	
							maxCloseScore = jumpScore;
						}
						break;
					case J_FAR:
						if (jumpScore > maxFarScore)
						{
							extendsFar = true;
							maxFarId = pathId;
							maxFarScore = jumpScore;
						}
						break;
				}
			}
			//update the best close extension
			if (extendsClose)
			{
				eraseMarks.erase(maxCloseId);
				extPaths[maxCloseId].ovlp.curEnd = curPos;
				extPaths[maxCloseId].ovlp.extEnd = extPos;
				extPaths[maxCloseId].ovlp.score = maxCloseScore;
				extPaths[maxCloseId].shifts.push_back(curPos - extPos);

				if (_keepAlignment)
				{
					auto& kmerMatches = extPaths[maxCloseId].ovlp.kmerMatches;
					if (kmerMatches.empty() || curPos - kmerMatches.back().first > 
											   (int32_t)Parameters::get().kmerSize)
					{
						kmerMatches.emplace_back(curPos, extPos);
					}
				}
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().ovlp.curEnd = curPos;
				extPaths.back().ovlp.extEnd = extPos;
				extPaths.back().ovlp.score = maxFarScore;
				extPaths.back().shifts.push_back(curPos - extPos);

				if (_keepAlignment)
				{
					auto& kmerMatches = extPaths.back().ovlp.kmerMatches;
					if (kmerMatches.empty() || curPos - kmerMatches.back().first > 
											   (int32_t)Parameters::get().kmerSize)
					{
						kmerMatches.emplace_back(curPos, extPos);
					}
				}
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				(this->goodStart(curPos, extPos, curLen, extLen) ||
				fastaRec.id == extReadPos.readId.rc()))	//FIXME: temporary bypass overhang
			{
				OverlapRange ovlp(fastaRec.id, extReadPos.readId,
								  curPos, extPos, curLen, extLen);
				extPaths.emplace_back(ovlp);
			}
			//cleaning up
			for (auto itEraseId = eraseMarks.rbegin(); 
				 itEraseId != eraseMarks.rend(); ++itEraseId)
			{
				extPaths[*itEraseId] = extPaths.back();
				extPaths.pop_back();
			}
		} //end loop over kmer occurences in other reads
	} //end loop over kmers in the current read

	//copy to coplete paths
	for (auto& ap : activePaths)
	{
		for (auto& dpRec : ap.second)
		{
			if (this->overlapTest(dpRec.ovlp))
			{
				completedPaths[ap.first].push_back(dpRec);
			}
		}
	}

	//leave only one overlap for each starting position
	for (auto& ap : completedPaths)
	{
		std::unordered_map<std::pair<int32_t, int32_t>, 
						   DPRecord, pairhash> maxByStart;
		for (auto& dpRec : ap.second)
		{
			auto& curMax = maxByStart[std::make_pair(dpRec.ovlp.curBegin, 
													 dpRec.ovlp.extBegin)];
			if (dpRec.ovlp.score > curMax.ovlp.score)
			{
				curMax = dpRec;
			}
		}
		ap.second.clear();
		for (auto& maxDp : maxByStart) ap.second.push_back(maxDp.second);
	}

	//copmutes sequence divergence rate
	auto computeSequenceDiv = [this, &fastaRec](DPRecord& rec)
	{
		//kmer divergence
		/*std::unordered_set<Kmer> curKmers;
		std::unordered_set<Kmer> extKmers;
		std::unordered_set<Kmer> allKmers;
		std::unordered_set<Kmer> sharedKmers;
		for (auto curKmerPos : IterKmers(fastaRec.sequence, rec.ovlp.curBegin,
										 rec.ovlp.curRange()))
		{
			curKmers.insert(curKmerPos.kmer);
			allKmers.insert(curKmerPos.kmer);
		}
		for (auto extKmerPos : IterKmers(_seqContainer.getSeq(rec.ovlp.extId),
										 rec.ovlp.extBegin, rec.ovlp.extRange()))
		{
			if (curKmers.count(extKmerPos.kmer)) sharedKmers.insert(extKmerPos.kmer);
			extKmers.insert(extKmerPos.kmer);
			allKmers.insert(extKmerPos.kmer);
		}
		float jaccardIndex = (float)sharedKmers.size() / allKmers.size();
		float kmerDiv = -1.0f / Parameters::get().kmerSize * 
						log(2 * jaccardIndex / (1 + jaccardIndex));*/
		//float kmerDiv = ((float)sharedKmers.size() / curKmers.size() + 
		//				 (float)sharedKmers.size() / extKmers.size()) / 2;
		//float kmerDistance = 1 - pow(kmerDiv, 1.0 / Parameters::get().kmerSize);

		std::vector<uint8_t> dnaCur(rec.ovlp.curRange());
		for (size_t i = 0; i < (size_t)rec.ovlp.curRange(); ++i)
		{
			dnaCur[i] = fastaRec.sequence.atRaw(rec.ovlp.curBegin + i);
		}
		std::vector<uint8_t> dnaExt(rec.ovlp.extRange());
		auto& extSeq = _seqContainer.getSeq(rec.ovlp.extId);
		for (size_t i = 0; i < (size_t)rec.ovlp.extRange(); ++i)
		{
			dnaExt[i] = extSeq.atRaw(rec.ovlp.extBegin + i);
		}

		int seqDiff = abs((int)dnaCur.size() - (int)dnaExt.size());
		int bandWidth = seqDiff + 100;
		float kswDiv = kswAlign(dnaCur, dnaExt, 1, -2, 2, 1, bandWidth);
		rec.ovlp.divergence = kswDiv;

		//std::cout << rec.ovlp.curRange()
		//	<< "\t" << div1500 << "\t" << div500 << "\t" << kswDiv << std::endl;

	};
	
	std::vector<OverlapRange> detectedOverlaps;
	auto addOverlap = [curLen, computeSequenceDiv, &detectedOverlaps, this]
		(DPRecord& rec, int32_t extLen)
	{
		rec.ovlp.leftShift = median(rec.shifts);
		rec.ovlp.rightShift = extLen - curLen + 
								rec.ovlp.leftShift;
		computeSequenceDiv(rec);
		if (rec.ovlp.divergence < _overlapDivergenceThreshold)
		{
			detectedOverlaps.push_back(rec.ovlp);
		}
	};

	for (auto& ap : completedPaths)
	{
		int32_t extLen = _seqContainer.seqLen(ap.first);
		DPRecord* maxRecord = nullptr;
		bool passedTest = false;
		for (auto& dpRec : ap.second)
		{
			if (!uniqueExtensions)
			{
				addOverlap(dpRec, extLen);
			}
			else
			{
				passedTest = true;
				if (!maxRecord || dpRec.ovlp.score > maxRecord->ovlp.score)
				{
					maxRecord = &dpRec;
				}
			}
		}
		if (uniqueExtensions && passedTest)
		{
			addOverlap(*maxRecord, extLen);
		}
	}

	for (auto& ovlp : detectedOverlaps)
	{
		if (_dumpFile) _dumpFile << ovlp.divergence << "\n";
	}

	return detectedOverlaps;
}


float OverlapContainer::meanDivergence()
{
	double sumDiv = 0.0f;
	int numDiv = 0;
	for (auto& seqOvlps : _overlapIndex)
	{
		for (auto& ovlp : seqOvlps.second)
		{
			sumDiv += ovlp.divergence;
			++numDiv;
		}
	}
	return numDiv ? sumDiv / numDiv : 0.0f;
}

std::vector<OverlapRange>
OverlapContainer::seqOverlaps(FastaRecord::Id seqId) const
{
	const FastaRecord& record = _queryContainer.getIndex().at(seqId);
	return _ovlpDetect.getSeqOverlaps(record, _onlyMax);
}

std::vector<OverlapRange>
	OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	_indexMutex.lock();
	if (!_cached.count(readId))
	{
		_indexMutex.unlock();
		auto overlaps = this->seqOverlaps(readId);
		_indexMutex.lock();
		this->storeOverlaps(overlaps, readId);
	}
	auto overlaps = _overlapIndex.at(readId);
	_indexMutex.unlock();
	return overlaps;
}

void OverlapContainer::storeOverlaps(const std::vector<OverlapRange>& overlaps,
									 FastaRecord::Id seqId)
{
	_cached.insert(seqId);
	_cached.insert(seqId.rc());

	auto& fwdOverlaps = _overlapIndex[seqId];
	auto& revOverlaps = _overlapIndex[seqId.rc()];

	std::unordered_set<FastaRecord::Id> extisting;
	if (_onlyMax)
	{
		for (auto& ovlp : fwdOverlaps) extisting.insert(ovlp.extId);
	}

	for (auto& ovlp : overlaps)
	{
		if (_onlyMax && extisting.count(ovlp.extId)) continue;

		auto revOvlp = ovlp.reverse();
		fwdOverlaps.push_back(ovlp);
		revOverlaps.push_back(ovlp.complement());
		_overlapIndex[revOvlp.curId].push_back(revOvlp);
		_overlapIndex[revOvlp.curId.rc()].push_back(revOvlp.complement());
	}
}

void OverlapContainer::findAllOverlaps()
{
	//Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (auto& hashPair : _queryContainer.getIndex())
	{
		allQueries.push_back(hashPair.first);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, &indexMutex] (const FastaRecord::Id& seqId)
	{
		auto& fastaRec = _queryContainer.getIndex().at(seqId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(fastaRec, false);

		indexMutex.lock();
		this->storeOverlaps(overlaps, seqId);
		indexMutex.unlock();
	};

	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);

	int numOverlaps = 0;
	for (auto& seqOvlps : _overlapIndex) numOverlaps += seqOvlps.second.size();
	Logger::get().debug() << "Found " << numOverlaps << " overlaps";

	this->filterOverlaps();


	numOverlaps = 0;
	for (auto& seqOvlps : _overlapIndex) numOverlaps += seqOvlps.second.size();
	Logger::get().debug() << "Left " << numOverlaps 
		<< " overlaps after filtering";
}


void OverlapContainer::filterOverlaps()
{
	std::vector<FastaRecord::Id> seqIds;
	for (auto seqIndex : _queryContainer.getIndex())
	{
		seqIds.push_back(seqIndex.first);
	}

	std::function<void(const FastaRecord::Id& seqId)> filterParallel =
	[this] (const FastaRecord::Id& seqId)
	{
		auto& overlaps = _overlapIndex[seqId];
		
		std::vector<SetNode<OverlapRange*>*> overlapSets;
		for (auto& ovlp : overlaps) 
		{
			overlapSets.push_back(new SetNode<OverlapRange*>(&ovlp));
		}
		for (size_t i = 0; i < overlapSets.size(); ++i)
		{
			for (size_t j = 0; j < overlapSets.size(); ++j)
			{
				OverlapRange& ovlpOne = *overlapSets[i]->data;
				OverlapRange& ovlpTwo = *overlapSets[j]->data;

				if (ovlpOne.extId != ovlpTwo.extId) continue;
				int curDiff = ovlpOne.curRange() - ovlpOne.curIntersect(ovlpTwo);
				int extDiff = ovlpOne.extRange() - ovlpOne.extIntersect(ovlpTwo);

				if (curDiff < 100 && extDiff < 100) 
				{
					unionSet(overlapSets[i], overlapSets[j]);
				}
			}
		}
		std::unordered_map<SetNode<OverlapRange*>*, 
						   std::vector<OverlapRange>> clusters;
		for (auto& ovlp: overlapSets) 
		{
			clusters[findSet(ovlp)].push_back(*ovlp->data);
		}
		overlaps.clear();
		for (auto& cluster : clusters)
		{
			OverlapRange* maxOvlp = nullptr;
			for (auto& ovlp : cluster.second)
			{
				if (!maxOvlp || ovlp.score > maxOvlp->score)
				{
					maxOvlp = &ovlp;
				}
			}
			overlaps.push_back(*maxOvlp);
		}
		for (auto& ovlpNode : overlapSets) delete ovlpNode;

		std::sort(overlaps.begin(), overlaps.end(), 
				  [](const OverlapRange& o1, const OverlapRange& o2)
				  {return o1.curBegin < o2.curBegin;});

	};
	processInParallel(seqIds, filterParallel, 
					  Parameters::get().numThreads, false);
}


/*void OverlapContainer::saveOverlaps(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	for (auto& hashPair : _overlapIndex)
	{
		for (auto& ovlp : hashPair.second)
		{
			fout << ovlp.serialize() << std::endl;
		}
	}
}

void OverlapContainer::loadOverlaps(const std::string& filename)
{
	std::ifstream fin(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	std::string buffer;
	while(!fin.eof())
	{
		std::getline(fin, buffer);
		if (buffer.empty()) break;
		OverlapRange ovlp;
		ovlp.unserialize(buffer);
		_overlapIndex[ovlp.curId].push_back(ovlp);
	}
}*/

void OverlapContainer::buildIntervalTree()
{
	//Logger::get().debug() << "Building interval tree";
	for (auto& seqOvlps : _overlapIndex)
	{
		std::vector<Interval<OverlapRange*>> intervals;
		for (auto& ovlp : seqOvlps.second)
		{
			intervals.emplace_back(ovlp.curBegin, ovlp.curEnd, &ovlp);
		}
		_ovlpTree[seqOvlps.first] = IntervalTree<OverlapRange*>(intervals);
	}
}

std::vector<Interval<OverlapRange*>> 
	OverlapContainer::getOverlaps(FastaRecord::Id seqId, 
								  int32_t start, int32_t end) const
{
	return _ovlpTree.at(seqId).findOverlapping(start, end);
}

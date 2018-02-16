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
#include "../common/matrix.h"

namespace
{
	//banded glocal alignment
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   std::string& outOne, std::string& outTwo, 
						   int bandWidth)
	{
		static const int32_t MATCH = 5;
		static const int32_t SUBST = -5;
		static const int32_t INDEL = -3;

		static const int32_t matchScore[] = {SUBST, MATCH};
		static const int32_t indelScore[] = {INDEL, 0};

		const size_t numRows = seqOne.length() + 1;
		const size_t numCols = 2 * bandWidth + 1;

		Matrix<char> backtrackMatrix(numRows, numCols);
		std::vector<int32_t> scoreRowOne(numCols, 0);
		std::vector<int32_t> scoreRowTwo(numCols, 0);


		for (size_t i = 0; i < numRows; ++i) 
		{
			size_t j = std::max(0, bandWidth - (int)i);
			backtrackMatrix.at(i, j) = 1;
		}
		for (size_t i = 0; i < numCols; ++i) 
		{
			backtrackMatrix.at(0, i) = 0;
		}

		//filling DP matrices
		for (size_t i = 1; i < numRows; ++i)
		{
			int leftOverhang = bandWidth - (int)i + 1;
			int rightOverhand = (int)i + bandWidth - (int)seqTwo.length();
			size_t colLeft = std::max(0, leftOverhang);
			size_t colRight = std::min((int)numCols, (int)numCols - rightOverhand);

			for (int j = colLeft; j < (int)colRight; ++j)
			{
				size_t twoCoord = j + i - bandWidth;
				int32_t cross = scoreRowOne[j] + 
								matchScore[seqOne[i - 1] == seqTwo[twoCoord - 1]];
				char maxStep = 2;
				int32_t maxScore = cross;

				if (j < (int)numCols - 1) //up
				{
					int32_t up = scoreRowOne[j + 1] +
								 indelScore[twoCoord == seqTwo.length()];
					if (up > maxScore)
					{
						maxStep = 1;
						maxScore = up;
					}
				}

				if (j > 0) //left
				{
					int32_t left = scoreRowTwo[j - 1] + 
								   indelScore[i == seqOne.length()];
					if (left > maxScore)
					{
						maxStep = 0;
						maxScore = left;
					}
				}

				scoreRowTwo[j] = maxScore;
				backtrackMatrix.at(i, j) = maxStep;
			}
			scoreRowOne.swap(scoreRowTwo);
		}

		//backtrack
		outOne.reserve(seqOne.length() * 3 / 2);
		outTwo.reserve(seqTwo.length() * 3 / 2);

		int i = numRows - 1;
		int j = bandWidth - (int)numRows + (int)seqTwo.length() + 1;
		while (i != 0 || j != bandWidth) 
		{
			size_t twoCoord = j + i - bandWidth;
			if(backtrackMatrix.at(i, j) == 1) //up
			{
				outOne += seqOne[i - 1];
				outTwo += '-';
				i -= 1;
				j += 1;
			}
			else if (backtrackMatrix.at(i, j) == 0) //left
			{
				outOne += '-';
				outTwo += seqTwo[twoCoord - 1];
				j -= 1;
			}
			else	//cross
			{
				outOne += seqOne[i - 1];
				outTwo += seqTwo[twoCoord - 1];
				i -= 1;
			}
		}
		std::reverse(outOne.begin(), outOne.end());
		std::reverse(outTwo.begin(), outTwo.end());
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
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
				< _maxJump / CLOSE_JUMP &&
			curNext - curPrev < _gapSize)
		{
			return J_CLOSE;
		}

		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maxJump / FAR_JUMP)
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

	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength / OVLP_DIVERGENCE)
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
		//test with rela pairwise alignment
		std::string seqA = fastaRec.sequence.substr(rec.ovlp.curBegin, 
													rec.ovlp.curRange()).str();
		std::string seqB = _seqContainer.getSeq(rec.ovlp.extId)
								.substr(rec.ovlp.extBegin, 
										rec.ovlp.extRange()).str();
		std::string alnA;
		std::string alnB;
		const int bandWidth = abs((int)seqA.length() - 
								  (int)seqB.length()) + 
								  		Config::get("maximum_jump");
		//std::cout << ">>>>>" << std::endl;
		pairwiseAlignment(seqA, seqB, alnA, alnB, bandWidth);
		int matches = 0;
		for (size_t i = 0; i < alnA.size(); ++i)
		{
			if (alnA[i] == alnB[i]) ++matches;
		}
		float trueDiv = 1 - (float)matches/ alnA.size();

		//kmer divergence
		std::unordered_set<Kmer> curKmers;
		std::unordered_set<Kmer> extKmers;
		std::unordered_set<Kmer> sharedKmers;
		for (auto curKmerPos : IterKmers(fastaRec.sequence, rec.ovlp.curBegin,
										 rec.ovlp.curRange()))
		{
			curKmers.insert(curKmerPos.kmer);
		}
		for (auto extKmerPos : IterKmers(_seqContainer.getSeq(rec.ovlp.extId),
										 rec.ovlp.extBegin, rec.ovlp.extRange()))
		{
			if (curKmers.count(extKmerPos.kmer)) sharedKmers.insert(extKmerPos.kmer);
			extKmers.insert(extKmerPos.kmer);
		}

		float kmerDiv = ((float)sharedKmers.size() / curKmers.size() + 
						 (float)sharedKmers.size() / extKmers.size()) / 2;
		float kmerDistance = 1 - pow(kmerDiv, 1.0 / Parameters::get().kmerSize);
		
		if (fabsf(trueDiv - kmerDistance) > 0.05)
		{
			for (size_t i = 0; i < alnA.size() / 100 + 1; ++i)
			{
				std::cout << alnA.substr(i * 100, 100) << "\n" 
					<< alnB.substr(i * 100, 100) << "\n\n";
			}
		}
		std::cout << rec.ovlp.curRange() << "\t" << trueDiv 
			<< "\t" << kmerDistance << std::endl;
		rec.ovlp.divergence = kmerDistance;
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



//pre-computes sequence-to-kmer divergence function by
//direct simulation
/*void OverlapDetector::preComputeKmerDivergence()
{
	const int NUM_BASES = 1000000;
	auto simulate = [](float seqDivergence, int numBases)
	{
		std::vector<char> bases(numBases, true);
		for (size_t i = 0; i < (size_t)(seqDivergence * numBases); ++i)
		{
			size_t pos = random() % (NUM_BASES - Parameters::get().kmerSize);
			for (size_t j = 0; j < Parameters::get().kmerSize; ++j)
			{
				bases[pos + j] = 0;
			}
		}
		size_t corrupted = false;
		for (auto c : bases) if (!c) ++corrupted;
		return (float)corrupted / numBases;
	};

	for (size_t i = 0; i < 50; ++i)
	{
		float seqDiv = (float)i / 100;
		float kmerDiv = simulate(seqDiv, NUM_BASES);
		_seqToKmerDiv.push_back(kmerDiv);
	}
}

float OverlapDetector::kmerToSeqDivergence(float kmerDivergence) const
{
	float lowestDiff = std::numeric_limits<float>::max();
	float seqDiv = 0.0f;
	for (size_t i = 0; i < _seqToKmerDiv.size(); ++i)
	{
		if (fabs(_seqToKmerDiv[i] - kmerDivergence) < lowestDiff)
		{
			lowestDiff = fabs(_seqToKmerDiv[i] - kmerDivergence);
			seqDiv = (float)i / 100;
		}
	}
	return seqDiv;
}*/

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

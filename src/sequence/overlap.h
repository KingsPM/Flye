//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unordered_set>
#include <mutex>
#include <sstream>
#include <fstream>

#include <cuckoohash_map.hh>
#include "IntervalTree.h"

#include "vertex_index.h"
#include "sequence_container.h"
#include "../common/logger.h"
#include "../common/progress_bar.h"


struct OverlapRange
{
	OverlapRange(FastaRecord::Id curId = FastaRecord::ID_NONE, 
				 FastaRecord::Id extId = FastaRecord::ID_NONE, 
				 int32_t curInit = 0, int32_t extInit = 0,
				 int32_t curLen = 0, int32_t extLen = 0): 
		curId(curId), curBegin(curInit), curEnd(curInit), curLen(curLen),
		extId(extId), extBegin(extInit), extEnd(extInit), extLen(extLen),
		leftShift(0), rightShift(0), score(0), divergence(0.0f)
	{}
	int32_t curRange() const {return curEnd - curBegin;}
	int32_t extRange() const {return extEnd - extBegin;}

	OverlapRange reverse() const
	{
		OverlapRange rev(*this);
		std::swap(rev.curId, rev.extId);
		std::swap(rev.curBegin, rev.extBegin);
		std::swap(rev.curEnd, rev.extEnd);
		std::swap(rev.curLen, rev.extLen);
		rev.leftShift = -rev.leftShift;
		rev.rightShift = -rev.rightShift;

		for (auto& posPair : rev.kmerMatches) 
		{
			std::swap(posPair.first, posPair.second);
		}
		std::sort(rev.kmerMatches.begin(), rev.kmerMatches.end(),
				  [](const std::pair<int32_t, int32_t>& p1,
				  	 const std::pair<int32_t, int32_t>& p2)
				  	 {return p1.first < p2.first;});

		return rev;
	}

	OverlapRange complement() const
	{
		OverlapRange comp(*this);
		std::swap(comp.leftShift, comp.rightShift);
		comp.leftShift = -comp.leftShift;
		comp.rightShift = -comp.rightShift;

		std::swap(comp.curBegin, comp.curEnd);
		comp.curBegin = curLen - comp.curBegin - 1;
		comp.curEnd = curLen - comp.curEnd - 1;

		std::swap(comp.extBegin, comp.extEnd);
		comp.extBegin = extLen - comp.extBegin - 1;
		comp.extEnd = extLen - comp.extEnd - 1;

		comp.curId = comp.curId.rc();
		comp.extId = comp.extId.rc();

		for (auto& posPair : comp.kmerMatches) 
		{
			posPair.first = curLen - posPair.first - 1, 
			posPair.second = extLen - posPair.second - 1;
		}
		std::reverse(comp.kmerMatches.begin(), comp.kmerMatches.end());

		return comp;
	}

	int32_t project(int32_t curPos) const
	{
		if (kmerMatches.empty())
		{
			float lengthRatio = (float)this->extRange() / this->curRange();
			int32_t projectedPos = extBegin + 
							float(curPos - curBegin) * lengthRatio;
			return std::max(extBegin, std::min(projectedPos, extEnd));
		}
		else
		{
			auto cmpFirst = [] (const std::pair<int32_t, int32_t>& pair, 
							  	int32_t value)
								{return pair.first < value;};
			auto iter = std::lower_bound(kmerMatches.begin(), kmerMatches.end(), 
										 curPos, cmpFirst);
			if (iter == kmerMatches.end()) return kmerMatches.back().second;

			//Logger::get().debug() << curPos - curBegin << " " << iter->second - extBegin <<
			//	" " << curPos << " " << curBegin << " " << curEnd << " " << kmerMatches.size();
			return iter->second;
		}
	}

	bool contains(int32_t curPos, int32_t extPos) const
	{
		return curBegin <= curPos && curPos <= curEnd &&
			   extBegin <= extPos && extPos <= extEnd;
	}

	bool containedBy(const OverlapRange& other) const
	{
		if (curId != other.curId || extId != other.curId) return false;

		return other.curBegin < curBegin && curEnd < other.curEnd &&
			   other.extBegin < extBegin && extEnd <= other.extEnd;
	}

	int32_t curIntersect(const OverlapRange& other) const
	{
		return std::min(curEnd, other.curEnd) - 
			   std::max(curBegin, other.curBegin);
	}

	int32_t extIntersect(const OverlapRange& other) const
	{
		return std::min(extEnd, other.extEnd) - 
			   std::max(extBegin, other.extBegin);
	}

	bool equals(const OverlapRange& other) const
	{
		return other.curId == curId && other.extId == extId &&
			   other.curBegin == curBegin && other.curEnd == curEnd &&
			   other.extBegin == extBegin && other.extEnd == extEnd;
	}

	std::string serialize() const
	{
		std::stringstream ss;
		ss << curId << " " << curBegin << " " << curEnd << " " 
		   << leftShift << " " << extId << " " << extBegin << " " 
		   << extEnd << " " << rightShift;
		return ss.str();
	}

	void unserialize(const std::string& str)
	{
		std::stringstream ss(str);
		ss >> curId >> curBegin >> curEnd >> leftShift 
		   >> extId >> extBegin >> extEnd >> rightShift;
	}

	//current read
	FastaRecord::Id curId;
	int32_t curBegin;
	int32_t curEnd;
	int32_t curLen;

	//extension read
	FastaRecord::Id extId;
	int32_t extBegin;
	int32_t extEnd;
	int32_t extLen;

	int32_t leftShift;
	int32_t rightShift;

	int32_t score;
	float divergence;
	std::vector<std::pair<int32_t, int32_t>> kmerMatches;
};

class OverlapDetector
{
public:
	OverlapDetector(const SequenceContainer& seqContainer,
					const VertexIndex& vertexIndex,
					int maxJump, int minOverlap, int maxOverhang,
					int gapSize, float overlapDivThreshold, bool keepAlignment):
		_maxJump(maxJump),
		_minOverlap(minOverlap),
		_maxOverhang(maxOverhang),
		_gapSize(gapSize),
		_overlapDivergenceThreshold(overlapDivThreshold),
		_checkOverhang(maxOverhang > 0),
		_keepAlignment(keepAlignment),
		_vertexIndex(vertexIndex),
		_seqContainer(seqContainer)
	{}

	std::vector<OverlapRange> 
	getSeqOverlaps(const FastaRecord& fastaRec, bool uniqueExtensions) const;
	void setDumpFile(const std::string& file)
	{
		_dumpFile.open(file);
	}

private:
	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};

	
	bool    goodStart(int32_t currentPos, int32_t extensionPos, 
				      int32_t curLen, int32_t extLen) const;
	bool    overlapTest(const OverlapRange& ovlp) const;
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext) const;

	mutable std::ofstream _dumpFile;
	const int _maxJump;
	const int _minOverlap;
	const int _maxOverhang;
	const int _gapSize;
	const float _overlapDivergenceThreshold;
	const bool _checkOverhang;
	const bool _keepAlignment;

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;
};


class OverlapContainer
{
public:
	OverlapContainer(const OverlapDetector& ovlpDetect,
					 const SequenceContainer& queryContainer,
					 bool onlyMax):
		_ovlpDetect(ovlpDetect),
		_queryContainer(queryContainer),
		_onlyMax(onlyMax)
	{}

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;

	//void saveOverlaps(const std::string& filename);
	//void loadOverlaps(const std::string& filename);
	float meanDivergence();

	void findAllOverlaps();
	std::vector<OverlapRange> seqOverlaps(FastaRecord::Id readId) const;
	std::vector<OverlapRange> lazySeqOverlaps(FastaRecord::Id readId);
	const OverlapIndex& getOverlapIndex() const {return _overlapIndex;}

	void buildIntervalTree();
	std::vector<Interval<OverlapRange*>> 
		getOverlaps(FastaRecord::Id seqId, int32_t start, int32_t end) const;

private:
	void storeOverlaps(const std::vector<OverlapRange>& overlaps, 
					   FastaRecord::Id seqId);
	void filterOverlaps();

	const OverlapDetector& _ovlpDetect;
	const SequenceContainer& _queryContainer;
	const bool _onlyMax;

	std::mutex _indexMutex;
	OverlapIndex _overlapIndex;
	std::unordered_set<FastaRecord::Id> _cached;

	std::unordered_map<FastaRecord::Id, 
					   IntervalTree<OverlapRange*>> _ovlpTree;
};

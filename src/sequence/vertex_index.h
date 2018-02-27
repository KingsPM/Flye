//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <iostream>

#include <cuckoohash_map.hh>

#include "kmer.h"
#include "sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"

class VertexIndex
{
public:
	~VertexIndex()
	{
		this->clear();
	}
	VertexIndex(const SequenceContainer& seqContainer, int sampleRate):
		_seqContainer(seqContainer), _outputProgress(false), _sampleRate(sampleRate)
	{}

	VertexIndex(const VertexIndex&) = delete;
	void operator=(const VertexIndex&) = delete;

private:
	struct ReadPosition
	{
		ReadPosition(FastaRecord::Id readId = FastaRecord::ID_NONE, 
					 int32_t position = 0):
			readId(readId), position(position) {}
		FastaRecord::Id readId;
		int32_t position;
	};

	struct ReadVector
	{
		uint32_t capacity;
		uint32_t size;
		ReadPosition* data;
	};

public:
	typedef std::map<size_t, size_t> KmerDistribution;

	class KmerPosIterator
	{
	public:
		KmerPosIterator(ReadVector rv, size_t index, bool revComp, 
						const SequenceContainer& seqContainer):
			rv(rv), index(index), revComp(revComp), 
			seqContainer(seqContainer) 
		{}

		bool operator==(const KmerPosIterator& other) const
		{
			return rv.data == other.rv.data && index == other.index;
		}
		bool operator!=(const KmerPosIterator& other) const
		{
			return !(*this == other);
		}

		ReadPosition operator*() const
		{
			ReadPosition pos = rv.data[index];
			if (revComp)
			{
				pos.readId = pos.readId.rc();
				int32_t seqLen = seqContainer.seqLen(pos.readId);
				pos.position = seqLen - pos.position - 
							   Parameters::get().kmerSize;
			}
			return pos;
		}

		KmerPosIterator& operator++()
		{
			++index;
			return *this;
		}
	
	private:
		ReadVector rv;
		size_t index;
		bool revComp;
		const SequenceContainer& seqContainer;
	};

	class IterHelper
	{
	public:
		IterHelper(ReadVector rv, bool revComp, 
				   const SequenceContainer& seqContainer): 
			rv(rv), revComp(revComp), seqContainer(seqContainer) {}

		KmerPosIterator begin()
		{
			return KmerPosIterator(rv, 0, revComp, seqContainer);
		}

		KmerPosIterator end()
		{
			return KmerPosIterator(rv, rv.size, revComp, seqContainer);
		}

	private:
		ReadVector rv;
		bool revComp;
		const SequenceContainer& seqContainer;
	};


	void countKmers(size_t hardThreshold, int genomeSize);
	void buildIndex(int minCoverage, int maxCoverage);
	void clear();

	IterHelper iterKmerPos(Kmer kmer) const
	{
		bool revComp = kmer.standardForm();
		return IterHelper(_kmerIndex[kmer], revComp,
						  _seqContainer);
	}

	bool isSolid(Kmer kmer) const
	{
		kmer.standardForm();
		return _kmerIndex.contains(kmer);
	}

	void outputProgress(bool set) 
	{
		_outputProgress = set;
	}

	const KmerDistribution& getKmerHist() const
	{
		return _kmerDistribution;
	}

	int getSampleRate() const {return _sampleRate;}

private:
	void addFastaSequence(const FastaRecord& fastaRecord);

	const SequenceContainer& _seqContainer;
	KmerDistribution _kmerDistribution;
	bool _outputProgress;
	int32_t _sampleRate;

	const size_t INDEX_CHUNK = 32 * 1024 * 1024 / sizeof(ReadPosition);
	std::vector<ReadPosition*> _memoryChunks;
	cuckoohash_map<Kmer, ReadVector> _kmerIndex;
	cuckoohash_map<Kmer, size_t> _kmerCounts;
};

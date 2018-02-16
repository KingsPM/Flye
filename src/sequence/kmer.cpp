//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <cassert>

#include "kmer.h"
#include "vertex_index.h"
#include "sequence_container.h"
#include "../common/config.h"


Kmer::Kmer(const DnaSequence& dnaString, 
		   size_t start, size_t length):
	_representation(0)
{
	if (length != Parameters::get().kmerSize)
	{
		throw std::runtime_error("Kmer length inconsistency");
	}

	for (size_t i = start; i < start + length; ++i)	
	{
		_representation <<= 2;
		_representation += dnaString.atRaw(i);
	}
}

void Kmer::appendRight(DnaSequence::NuclType dnaSymbol)
{
	_representation <<= 2;
	_representation += dnaSymbol;

	KmerRepr kmerSize = Parameters::get().kmerSize;
	KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
	_representation &= kmerMask;
}

void Kmer::appendLeft(DnaSequence::NuclType dnaSymbol)
{
	_representation >>= 2;

	KmerRepr kmerSize = Parameters::get().kmerSize;
	KmerRepr shift = kmerSize * 2 - 2;
	_representation += dnaSymbol << shift;
}

Kmer Kmer::reverseComplement()
{
	KmerRepr tmpRepr = _representation;
	Kmer newKmer;

	for (unsigned int i = 0; i < Parameters::get().kmerSize; ++i)
	{
		newKmer._representation <<= 2;
		newKmer._representation += ~tmpRepr & 3;
		tmpRepr >>= 2;
	}

	return newKmer;
}

bool Kmer::standardForm()
{
	Kmer complKmer = this->reverseComplement();
	if (complKmer._representation < _representation)
	{
		_representation = complKmer._representation;
		return true;
	}
	return false;
}

KmerIterator::KmerIterator(const DnaSequence* readSeq, 
						   size_t position):
	_readSeq(readSeq),
	_position(position)
{
	if (position != readSeq->length() - Parameters::get().kmerSize)
	{
		//_kmer = Kmer(readSeq->substr(0, Parameters::get().kmerSize));
		_kmer = Kmer(*readSeq, 0, Parameters::get().kmerSize);
	}
}

bool KmerIterator::operator==(const KmerIterator& other) const
{
	return _readSeq == other._readSeq && _position == other._position;
}

bool KmerIterator::operator!=(const KmerIterator& other) const
{
	return !(*this == other);
}

KmerIterator& KmerIterator::operator++()
{
	size_t appendPos = _position + Parameters::get().kmerSize;
	_kmer.appendRight(_readSeq->atRaw(appendPos));
	++_position;
	return *this;
}

KmerPosition KmerIterator::operator*() const
{
	return KmerPosition(_kmer, _position);
}


KmerIterator IterKmers::begin()
{
	if (_sequence.length() < Parameters::get().kmerSize + _start) 
		return this->end();

	return KmerIterator(&_sequence, _start);
}

KmerIterator IterKmers::end()
{
	size_t end = _length == std::string::npos ? 
					_sequence.length() : _length + _start;
	return KmerIterator(&_sequence, end - Parameters::get().kmerSize);
}

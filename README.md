Flye assembler (successor of ABruijn)
=====================================

### Version: 2.3.2

Flye is a de novo assembler for long and noisy reads, such as
those produced by PacBio and Oxford Nanopore Technologies.
The algorithm uses an A-Bruijn graph to find the overlaps between reads
and does not require them to be error-corrected. After the initial assembly, 
Flye performs an extra repeat classification and analysis step to improve the structural
accuracy of the resulting sequence. The package also includes a polisher
module, which produces the final assembly of high nucleotide-level quality.


New in version 2.3
------------------

* ABruijn 2.x branch has been renamed to Flye, highlighting many substantial algorithmic changes
* Stable version of the repeat analysis module
* New command-line syntax (fallback mode with the old syntax is available)
* New --subassemblies mode for generating consensus of multiple assemblies
* Improved preformance and reduced memory footprint (now scales to human genome)
* Corrected reads are now supported
* Extra output with information about the contigs (coverage, multiplicity, graph paths etc.)
* Gzipped Fasta/q support
* Multiple read files support


Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)


Assembly graph
--------------

The Flye algorithms are operating on the assembly (repeat) graph. The edges in this graph 
represent genomic sequences, and nodes simply serve
as junctions. The genoimc chromosomes traverse this graph (in an unknown way) 
so as each unique edge is covered exactly once. The genomic repeats that were not
resolved are collapsed into the corresponding edges in the graph
(therefore genome structure remain umbigious).


<p align="center">
  <img src="docs/graph_example.png" alt="Graph example"/>
</p>

An example of a final assembly graph of a bacterial genome is above.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, while unique edges are black. The clusters of adjacent repeats are shown with the 
same color. Note that each edge is represented in two copies: forward and
reverse complement (marked with +/- signs), therefore the entire genome is
represented in two copies as well. Sometimes (as in this example), forward and reverse-complement
components are clearly separated, but often they form a single connected component
(in case if the genome contain unresolved inverted repeats).

In this example, there are two unresolved repeats: (i) a red repeat of multiplicity two
and length 35k and (ii) a green repeat cluster of multiplicity three and length 34k - 36k.
As the repeats remained unresolved, there are no reads in the dataset that cover
those repeats in full.


Third-party
-----------

Flye package includes some third-party software:

* [libcuckoo](http://github.com/efficient/libcuckoo)
* [intervaltree](https://github.com/ekg/intervaltree)
* [lemon](http://lemon.cs.elte.hu/trac/lemon)
* [minimap2](https://github.com/lh3/minimap2)


License
-------

Flye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Flye was developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Original assembler code: Yu Lin
* Original polisher code: Jeffrey Yuan
* Repeat graph and current package support: Mikhail Kolmogorov


Publications
------------
Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using Repeat Graphs", bioRxiv, 2018

Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Mark Chaisson and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs", PNAS, 2016


Contacts
--------
Please report any problems directly to the github issue tracker.
Also, you can send feedback to fenderglass@gmail.com

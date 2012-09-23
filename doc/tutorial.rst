========
Tutorial
========

This tutorial demonstrates the usage of *gat* with
a range of simple examples. It requires access to 
the UCSC web browser.

Overlap of CpG islands with histone modification tracks
-------------------------------------------------------

Assume we are interested in the question if CpG islands
and histone H3 lysine 4 trimethylation co-occur. We will
test this using data from the UCSC genome browser.

Data download
=============

To obtain data for this tutorial, go to the `UCSC
genome browser`_ and select *human* (assembly: *Mar. 2006
(NCBI36/hg18)*) as the species to work with. Go to the
`Table Browser` and download the following tracks or
tables:

+------------------+---------------+----------------------------------------+----------------------------------------+
|*Filename*        |*Group*        |*Track*                                 |Table                                   |
+------------------+---------------+----------------------------------------+----------------------------------------+
|cpg.bed.gz        |Regulation     |CpG islands                             |cpgIslandExt                            |
+------------------+---------------+----------------------------------------+----------------------------------------+
|h3k4m3.bed.gz     |Regulation     |Broad Histone                           |wgEncodeBroadChipSeqPeaksGm12878H3k4me3 |
+------------------+---------------+----------------------------------------+----------------------------------------+
|h3k27ac.bed.gz    |Regulation     |Broad Histone                           |wgEncodeBroadChipSeqPeaksGm12878H3k27ac |
+------------------+---------------+----------------------------------------+----------------------------------------+
|h3k4me1.bed.gz    |Regulation     |Broad Histone                           |wgEncodeBroadChipSeqPeaksGm12878H3k4me1 |
+------------------+---------------+----------------------------------------+----------------------------------------+

Save these tracks as bed-formatted files. 

Alternatively, the data can be obtained by running the following statements from
the command line::

   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -B -N -e "SELECT chrom, chromStart, chromEnd, 'cpg' FROM hg18.cpgIslandExt" | gzip > cpg.bed.gz
   for x in H3k4me3 H3k4me1 H3k27ac; do 
      mysql --user=genome --host=genome-mysql.cse.ucsc.edu -B -N -e "SELECT chrom, chromStart, chromEnd, '${x}' FROM hg18.wgEncodeBroadChipSeqPeaksGm12878${x}" > ${x}.bed.gz;
   done

Finally, we need a bed-formatted table with chromosome size::

   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -B -N -e "SELECT chrom, 0, size, 'workspace' FROM hg18.chromInfo" | gzip > genome.bed.gz

Running gat
===========

To run gat, execute::

   gat-run.py --segments=cpg.bed.gz --annotations="H3*.bed.gz" --workspace=genome.bed.gz --num-samples=1000 > gat1.tsv

The above command computes the overlap of cpg islands with any of the tracks in the files :file:`H3*.bed.gz`
within the full genome workspace. Logging information and the final results are written to :file:`stdout`, which 
is redirected to the file :file:`gat1.tsv`.

The contents of the output file :file:`gat1.tsv`:

+-----+----------+--------+------------+------------+------------+----------+-------+----------+----------+
|track|annotation|observed|expected    |CI95low     |CI95high    |stddev    |fold   |pvalue    |qvalue    |
+-----+----------+--------+------------+------------+------------+----------+-------+----------+----------+
|cpg  |H3k4me1   |4329604 |1862754.7410|1794054.0000|1929305.0000|41070.4838|2.3243 |1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+-------+----------+----------+
|cpg  |H3k27ac   |8674660 |975847.4750 |926041.0000 |1028620.0000|31680.2085|8.8894 |1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+-------+----------+----------+
|cpg  |H3k4me3   |12026180|953389.2270 |901977.0000 |1002005.0000|30086.4712|12.6141|1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+-------+----------+----------+

The first column contains the column headers. The second column, in words, states that 4,329,604 nucleotides 
within segments in track ``cpg`` overlap with segments ``H3k4me1``. In contrast, an overlap of only 1,862,755
is expected with a confidence interval of 1,794,054 and 1,929,305 (standard deviation: 41,070). This corresponds
to a 2.3 fold enrichment and significant of a level P < 0.001, the maximum obtainable from only 1,000 
randomizations.

Fold enrichments for the tracks ``H3k27ac`` and ``H3k4me3`` are 8.9 and 12.6, respectively. From this initial
analysis, we conclude that CpG islands indeed overlap with various regions of histone modification. But
is this the full story? Please read on.

Using workspaces
================

It is well known, that CpG islands often occur in or near promotors. It is also well known that 
various histone modifications have an effect of expression. Thus, the overlap between CpG islands
and the histone modification segments might simply be due to their association with genes. Using
the full genome as a background for the analysis might thus over-estimate their association.
Instead we should restrict the analysis to only those genomic segments that contain genes.

Here, we use refseq genes:

   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -B -N -e "SELECT chrom, GREATEST(0, txStart - 10000), txEnd + 10000, 'gene_workspace' FROM hg18.refGene" | gzip > genes.bed.gz 

The transcripts are extended by subtracting and adding 10,000 to the transcript start and end, respectively. 
Note that we do not have to worry about overlapping segments, as GAT will normalize each workspace before
analysis, i.e., overlapping segments within a workspace will be merged.

We can re-run the previous analysis by adding the additional workspace at the command line:

   gat-run.py --segments=cpg.bed.gz --annotations="H3*.bed.gz" --workspace=genome.bed.gz --workspace=genes.bed.gz --num-samples=1000 > gat2.tsv

Supplied with several workspaces, GAT will first normalize each workspace individually. It will then intersect all workspaces and restrict 
analysis to only those genomic segments, that are the intersection
of all workspaces. Providing both the genomic workspace :file:`genome.bed.gz` and the genes workspace :file:`genes.bed.gz` 
makes sure that transcript ends will be truncated if they extend beyond the chromosome end.

+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|track|annotation|observed|expected    |CI95low     |CI95high    |stddev    |fold  |pvalue    |qvalue    |
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|cpg  |H3k4me1   |3825525 |2512635.5170|2436271.0000|2592861.0000|48501.3396|1.5225|1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|cpg  |H3k27ac   |8270507 |1419942.7210|1361635.0000|1483282.0000|36962.9403|5.8245|1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|cpg  |H3k4me3   |11198479|1404928.3570|1343322.0000|1467763.0000|37930.9547|7.9709|1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+

Note how fold enrichment values have dropped, while results still remain significant.

As we are looking at gene-regulation, we might be interested to test of there is significant overlap within promotors only.
Here, we define as promotor the genomic segments 5kb upstream of a transcription start site. First, we create a bed-file
with promotor regions for our gene set::

   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -B -N -e "SELECT chrom, if(strand = '+', GREATEST(0, txStart - 5000), txend), if( strand = '+', txstart, txend + 5000), 'promotor_workspace' FROM hg18.refGene" | gzip > promotors.bed.gz

Then, we re-run gat by using the promotor workspace instead of the gene workspace::

   gat-run.py --segments=cpg.bed.gz --annotations="H3*.bed.gz" --workspace=genome.bed.gz --workspace=promotors.bed.gz --log=log --num-samples=1000 > gat3.tsv

+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|track|annotation|observed|expected    |CI95low     |CI95high    |stddev    |fold  |pvalue    |qvalue    |
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|cpg  |H3k4me1   |1741578 |1811927.6850|1765949.0000|1857663.0000|28290.4683|0.9612|8.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|cpg  |H3k27ac   |4006120 |1388385.4500|1342291.0000|1431366.0000|27043.2336|2.8855|1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+
|cpg  |H3k4me3   |5207804 |1769131.6190|1719615.0000|1816314.0000|30026.3604|2.9437|1.0000e-03|1.0000e+00|
+-----+----------+--------+------------+------------+------------+----------+------+----------+----------+

Fold enrichment values have dropped again. Enrichment with cpg islands with H3k27ac and H3k4me3 tracks
have converged to a value of approximately 2.9-fold enrichment. However, we now see a depletion of
H3k4me1 signals in cpg islands.

cpg islands: Methylation of CpG sites within the promoters of genes can lead to their silencing
In contrast, the hypomethylation of CpG sites has been associated with the over-expression of oncogenes

h3k4me1 activation
h3k27ac 
h3k4me3 activation


test::
   
   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -B -N -e "SELECT chrom, GREATEST(0, chromStart - 150000), chromEnd + 150000, trait FROM hg18.gwasCatalog" | gzip > gwas.bed.gz




.. UCSC genome browser: http://genome.ucsc.edu/
.. Table Browser: http://genome.ucsc.edu/cgi-bin/hgTables?org=Human&db=hg18



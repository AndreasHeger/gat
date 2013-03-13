SHELL=/bin/bash -O expand_aliases -l

SCRIPTSDIR=/home/andreas/devel/cgat/
SCRIPTSDIR=/ifs/devel/andreas/cgat/
# downloaded tabular results from GREAT server
# on 7.8.2012
# Options:
#   hg19
#   basal+extension rule (5kb/2kb + 1Mb)
#   only binomial test
GREAT_REFERENCE=srf.great.all.tsv

############################################################
# run a great analysis for functional annotation
great.tsv:
	gat-great.py --verbose=5 \
	        --log=$@.log \
	        --segments=srf.hg19.bed.gz \
		--annotations=regulatory_domains.bed.gz \
		--workspace=contigs_ungapped.bed.gz \
		--ignore-segment-tracks \
		--qvalue-method=BH \
		--descriptions=go2description.tsv \
	>& $@

############################################################
# run a gat analysis for functional annotation
# unnormed p-values
gatunnormed.tsv:
	gat-run.py --verbose=5 \
	        --log=$@.log \
	        --segments=srf.hg19.bed.gz \
		--annotations=regulatory_domains.bed.gz \
		--workspace=contigs_ungapped.bed.gz \
		--ignore-segment-tracks \
		--qvalue-method=BH \
		--descriptions=go2description.tsv \
	>& $@

############################################################
# run a gat analysis for functional annotation
# normed p-values
gatnormed.tsv:
	gat-run.py --verbose=5 \
	        --log=$@.log \
	        --segments=srf.hg19.bed.gz \
		--annotations=regulatory_domains.bed.gz \
		--workspace=contigs_ungapped.bed.gz \
		--ignore-segment-tracks \
		--qvalue-method=BH \
		--descriptions=go2description.tsv \
		--pvalue-method=norm \
	>& $@

############################################################
# run a gat analysis for functional annotation
gattruncated.tsv:
	gat-run.py --verbose=5 \
	        --log=$@.log \
	        --segments=srf.hg19.bed.gz \
		--annotations=regulatory_domains.bed.gz \
		--workspace=contigs_ungapped.bed.gz \
		--ignore-segment-tracks \
		--qvalue-method=BH \
		--descriptions=go2description.tsv \
		--pvalue-method=norm \
		--truncate-workspace-to-annotations \
	>& $@

############################################################
# run a gat analysis for functional annotation
# only two categories,
gatuntruncated_small.tsv:
	gat-run.py --verbose=5 \
	        --log=$@.log \
	        --segments=srf.hg19.bed.gz \
		--annotations=<( cat regulatory_domains.bed.gz | grep -e "GO:0030029" -e "GO:0015629") \
		--workspace=contigs_ungapped.bed.gz \
		--ignore-segment-tracks \
		--qvalue-method=BH \
		--descriptions=go2description.tsv \
		--pvalue-method=norm \
	>& $@

############################################################
# run a gat analysis for functional annotation
# only two categories, but with truncated workspace
gattruncated_small.tsv:
	gat-run.py --verbose=5 \
	        --log=$@.log \
	        --segments=srf.hg19.bed.gz \
		--annotations=<( cat regulatory_domains.bed.gz | grep -e "GO:0030029" -e "GO:0015629") \
		--workspace=contigs_ungapped.bed.gz \
		--ignore-segment-tracks \
		--qvalue-method=BH \
		--descriptions=go2description.tsv \
		--pvalue-method=norm \
		--truncate-workspace-to-annotations \
	>& $@

############################################################
greatserver.load: $(GREAT_REFERENCE)
	cat $< \
	| grep -e "# Ontology	" -e "^GO" \
	| perl -p -e "s/# //" \
	| csv2db.py --table=greatserver --index=ID > $@

############################################################
# ignore numpy warnings
%.load: %.tsv
	grep -v -e "__init__" -e "from cgat import" $< \
	| csv2db.py --table=$* --index=annotation > $@

############################################################
compare_great_fold.png:
	s3 csvdb "select ID, RegionFoldEnrich, fold from greatserver, great \
			WHERE great.annotation = greatserver.ID " \
	| python $(SCRIPTSDIR)plot_data.py --column=2,3 --no-legend \
		--log=$@.log \
		--xtitle="fold (GREAT)" \
		--ytitle="fold (GAT/GREAT)" \
		--logscale=xy \
		--hardcopy=$@

############################################################
compare_great_pvalue.png:
	s3 csvdb "select ID, BinomP, pvalue FROM greatserver, great \
			WHERE great.annotation = greatserver.ID " \
	| python $(SCRIPTSDIR)plot_data.py --column=2,3 --no-legend \
		--log=$@.log \
		--xtitle="pvalue (GREAT)" \
		--ytitle="pvalue (GAT/GREAT)" \
		--logscale=xy \
		--hardcopy=$@



compare_great_vs_gat_%_only_enrichment.png:
	s3 csvdb "SELECT a.$*, b.$* FROM gatunnormed as a, great as b \
			WHERE a.annotation = b.annotation AND a.fold > 1" \
	| python $(SCRIPTSDIR)plot_data.py \
		--hardcopy $@ \
		--logscale=xy \
		--xtitle="$* (GAT)" \
		--ytitle="$* (GREAT)"

compare_great_vs_gat_%_only_enrichment_truncated.png:
	s3 csvdb "SELECT a.$*, b.$* FROM gatunnormed as a, great as b \
			WHERE a.annotation = b.annotation AND a.fold > 1" \
	| python $(SCRIPTSDIR)plot_data.py \
		--hardcopy $@ \
		--logscale=xy \
		--xtitle="$* (GAT)" \
		--ytitle="$* (GREAT)" \
		--xrange=0.001,1 \
		--yrange=0.001,1

compare_great_vs_gat_%.png:
	s3 csvdb "SELECT a.$*, b.$* FROM gatunnormed as a, great as b WHERE a.annotation = b.annotation " \
	| python $(SCRIPTSDIR)plot_data.py \
		--hardcopy $@ \
		--logscale=xy \
		--xtitle="$* (GAT)" \
		--ytitle="$* (GREAT)"

compare_gattruncated_vs_gatnormed_%.png:
	s3 csvdb "SELECT a.$*, b.$* FROM gatnormed as a, gattruncated as b \
			WHERE a.annotation = b.annotation \
				AND a.pvalue not like '%nan%' AND b.pvalue not like '%nan%' " \
	| python $(SCRIPTSDIR)plot_data.py \
		--hardcopy $@ \
		--logscale=xy \
		--xtitle="$* (GAT normed)" \
		--ytitle="$* (GAT truncated)"

compare_great_vs_gatnormed_%_only_enrichment.png:
	s3 csvdb "SELECT a.$*, b.$* FROM gatnormed as a, great as b \
			WHERE a.annotation = b.annotation AND a.fold > 1" \
	| python $(SCRIPTSDIR)plot_data.py \
		--hardcopy $@ \
		--logscale=xy \
		--xtitle="$* (GAT)" \
		--ytitle="$* (GREAT)"


all: great.load gatunnormed.load greatserver.load gatnormed.load \
	gattruncated.load

png: compare_great_fold.png \
	compare_great_pvalue.png \
	compare_great_vs_gat_fold.png \
	compare_great_vs_gat_pvalue.png \
	compare_great_vs_gat_pvalue_only_enrichment.png \
	compare_great_vs_gat_pvalue_only_enrichment_truncated.png \
	compare_great_vs_gatnormed_pvalue_only_enrichment.png


clean:
	rm -f *.png *.load csvdb *.log gat*.tsv great.tsv 



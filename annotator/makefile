javac=javac

annotator:
	cd src; $(javac) -target 1.5 -source 1.5 -classpath ../lib/commons-cli-1.0.jar app/*.java segment/*.java
	cd src; jar -cvf ../lib/Annotator.jar app/*.class segment/*.class

distrib:
	tar -zchf Annotator.tgz makefile src/ doc/ lib/ scripts/ examples/*.sh examples/ann.ease.hg17.txt examples/id.hg17.txt examples/seg.hm.local.hg17.txt examples/syn.* examples/work.intron.cplusg.genscan.* examples/work.all*

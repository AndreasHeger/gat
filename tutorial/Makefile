all: TutorialFunctionalAnnotation.tar.gz TutorialGenomicAnnotation.tar.gz TutorialIntervalOverlap.tar.gz

TutorialFunctionalAnnotation.tar.gz:
	tar -cvzf $@ TutorialFunctionalAnnotation/*.bed.gz

TutorialGenomicAnnotation.tar.gz:
	tar -cvzf $@ TutorialGenomicAnnotation/{*.bed.gz,*.tsv}

TutorialIntervalOverlap.tar.gz:
	tar -cvzf $@ TutorialIntervalOverlap/*.bed.gz

publish:
	cp *.tar.gz /ifs/home/andreas/public_html/documentation/gat-examples/

.PHONY: clean article
article: ;
	@z 000_main&
	@z whatsleft&
clean:   ;
	@latexmk -silent -C > /dev/null
	@rm -f *.acn *.acr *.alg *.bbl *.glg *.glo *.gls *.ist *.run.xml

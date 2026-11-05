# Build the whole book, or a single chapter: `make`, `make ch06`, `make chapters`
LATEX = pdflatex -interaction=nonstopmode -halt-on-error
BUILD = build
CHAPTERS = $(patsubst chapters/%.tex,%,$(wildcard chapters/chapter*.tex))

.PHONY: book chapters clean
book:
	mkdir -p $(BUILD)
	$(LATEX) -output-directory=$(BUILD) structural_economics_book.tex
	$(LATEX) -output-directory=$(BUILD) structural_economics_book.tex

ch%: FORCE
	mkdir -p $(BUILD)
	cd chapters && $(LATEX) -output-directory=../$(BUILD) chapter$*.tex && $(LATEX) -output-directory=../$(BUILD) chapter$*.tex

chapters: $(CHAPTERS:chapter%=ch%)

clean:
	rm -rf $(BUILD)

FORCE:

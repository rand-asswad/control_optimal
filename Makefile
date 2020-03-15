doc: pdf html
	@echo "rendering document html and pdf..."

pdf: main.rmd
	@echo "rendering pdf book..."
	Rscript -e "bookdown::render_book('$<', 'bookdown::pdf_book', output_file='book.pdf')"
	[ -f _book/book.pdf ] && mv _book/book.pdf book.pdf
	@echo "done!"

html: main.rmd
	@echo "rendering html document..."
	Rscript -e "bookdown::render_book('$<', 'bookdown::html_document2')"
	@echo "done!"

clean:
	@echo "cleaning rmd output..."
	rm -rf _bookdown_files _book
slidify:
	R -e "slidify::slidify('index.Rmd')"

open:
	open index.html

clean:
	rm index.html index.md

clean_lib:
	rm -rf libraries

edit:
	subl index.Rmd

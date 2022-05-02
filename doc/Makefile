CPPFLAGS = -I../src/Engine -I../src/Primitives -I../src/Functions \
	   -I../../PsimagLite/src -I../../PsimagLite/src/Geometry 
CXX = g++

all: manual.pdf

FILES = $(shell find ../src -iname "*.h" -or -iname "*.cpp")
manual.tex: manual.ptex ../README.md  $(FILES)
	find ../src -iname "*.h" -or -iname "*.cpp" |\
              perl ../../PsimagLite/scripts/doc.pl manual.ptex

manual.pdf: manual.tex ../README.md  ../src/Engine/*h ../src/*cpp
	lualatex manual.tex
	bibtex manual.aux
	lualatex manual.tex
	lualatex manual.tex

clean:
	rm -f ../README.tex
	rm -f manual.tex manual.pdf



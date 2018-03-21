all: cmap
cmap: cmap.cc conto.cc ado.cc libfparser.h libfparser.cc
	c++ -o cmap cmap.cc libfparser.cc -lm
	strip cmap
clean:
	-rm cmap

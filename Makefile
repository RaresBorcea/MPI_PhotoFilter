.PHONY: build clean

build: filter

filter: filter.c
	mpicc -Wall filter.c -o filter
clean:
	rm -rf filter
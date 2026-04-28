CC ?= cc
CFLAGS ?= -O3 -std=c17
ARCH ?=
FLINT_CFLAGS := $(shell pkg-config --cflags flint)
FLINT_LIBS := $(shell pkg-config --libs flint)

.PHONY: all clean-c

all: hyperelliptic_finder_c

hyperelliptic_finder_c: hyperelliptic_finder.c
	$(CC) $(ARCH) $(CFLAGS) $< -o $@ $(FLINT_CFLAGS) $(FLINT_LIBS)

clean-c:
	rm -f hyperelliptic_finder_c

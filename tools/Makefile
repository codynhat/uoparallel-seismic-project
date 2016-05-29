#!/usr/bin/make -f

CC ?= cc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline
INCLUDE = -I../include
RM ?= rm -f

PROGS = velfileconvert

all: $(PROGS)

$(PROGS):%: %.c
	$(CC) $(CCFLAGS) $(INCLUDE) -o $@ $<

clean:
	$(RM) $(PROGS)

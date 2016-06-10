#!/bin/sh

VELFILE="../data/velocity-241-241-51.vbox"
STARTFILE="../docs/start-4.txt"
FORWARDSTARFILE="../docs/818-FS.txt"
OUTFILE="traveltime.ttbox"

./sweep-mpi-cuda $VELFILE $STARTFILE $FORWARDSTARFILE $OUTFILE

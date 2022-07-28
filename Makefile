# Makefile for HMMER
# 
##########
# HMMER - Biological sequence analysis with HMMs
# Copyright (C) 1992-1995 Sean R. Eddy
#
#   This source code is distributed under the terms of the 
#   GNU General Public License. See the files COPYING and 
#   LICENSE for details.
#    
###########

## where you want things installed.
#  In general BINDIR and MANDIR are all you need to think about 
#  changing in the Makefile.
#
BINDIR  = /usr/local/bin
MANDIR  = /usr/local/man
#BINDIR = $(HOME)/bin
#MANDIR = $(HOME)/man

#######
## You only need to change stuff below this if
## you experience compilation difficulties, or if you
## want to change the default compiler or compiler options.
#######

RELEASE     = "1.8.4"
RELEASEDATE = "July 1997"

## your compiler. This must be an ANSI-compliant compiler.
## In particular, SunOS cc is not ANSI. On SunOS, you should use gcc.
##
CC = cc
#CC = gcc	# GNU cc; especially on SunOS
#CC = fxc       # Alliant FX/2800 

## any special compiler flags you want
##
CFLAGS =  -O
#CFLAGS = -O -uniproc -w   # Alliant FX/2800

## machine specific definitions
# -DNOSTR                BSD?: no strstr() function
# -DNO_STRDUP            BSD?: no strdup() function
##
MDEFS =                    # Sun SPARC, Silicon Graphics IRIX, DEC Alpha
#MDEFS = -DNO_STRDUP       # Alliant FX/2800

## how to install the man pages 
## either cp -- to copy unformatted man page
## or a script with identical syntax to cp, to format & install formatted page
## Alternatively, you can set INSTMAN to something harmless like
## "echo", and hand-install the man pages.
##
INSTMAN   = cp             # this is fine on most machines 
MANSUFFIX = 1

#######
## You should not need to modify below this line.
#######
SHELL  = /bin/sh
LIBS   = -lm 

MANS  = hmmer hmma hmmb hmme hmmfs hmmls hmms hmmsw hmmt hmm_convert

PROGS = hmma hmmb hmme hmmfs hmmls hmms hmmsw hmmt hmm_convert
#       align_homologues is unsupported in this release
#       hmm_evolve is unsupported in this release

READMES = COPYING FILES GNULICENSE INSTALL README

HDRS =  config.h externs.h states.h version.h

SQUIDSRC = aligneval.c alignio.c sqerror.c sqio.c iupac.c msf.c revcomp.c\
	selex.c sre_ctype.c sre_math.c sre_string.c types.c cluster.c\
	weight.c stack.c dayhoff.c

SQUIDOBJ = aligneval.o alignio.o sqerror.o sqio.o iupac.o msf.o revcomp.o\
	selex.o sre_ctype.o sre_math.o sre_string.o types.o cluster.o\
	weight.o stack.o dayhoff.o

SRC =   align.c	build_main.c convert_main.c dbviterbi.c\
	align_homologues.c\
	emit.c forback.c fragviterbi.c hmma.c hmme.c hmmfs.c\
	hmmio.c	hmmls.c maxmodelmaker.c misc.c\
	output.c profiles.c prior.c\
	saviterbi.c scorestack.c search.c\
	states.c train_main.c viterbi.c\
	hmmsw.c	swviterbi.c trace.c weeviterbi.c

OBJ =   align.o dbviterbi.o emit.o forback.o fragviterbi.o hmmio.o \
	maxmodelmaker.o misc.o output.o prior.o profiles.o\
	saviterbi.o  scorestack.o states.o\
	viterbi.o swviterbi.o trace.o weeviterbi.o

SQUIDHDRS = squid.h sqfuncs.h


all: 	$(PROGS)

hmma: 	$(OBJ) $(SQUIDOBJ) hmma.o
	$(CC) $(CFLAGS) -o hmma hmma.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmmb:   $(OBJ) $(SQUIDOBJ) build_main.o
	$(CC) $(CFLAGS) -o hmmb  build_main.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmme:   $(OBJ) $(SQUIDOBJ) hmme.o
	$(CC) $(CFLAGS) -o hmme hmme.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmmfs:  $(OBJ) $(SQUIDOBJ) hmmfs.o
	$(CC) $(CFLAGS) -o hmmfs  hmmfs.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmmls:  $(OBJ) $(SQUIDOBJ) hmmls.o
	$(CC) $(CFLAGS) -o hmmls  hmmls.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmms:   $(OBJ) $(SQUIDOBJ) search.o
	$(CC) $(CFLAGS) -o hmms  search.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmmsw:  $(OBJ) $(SQUIDOBJ) hmmsw.o
	$(CC) $(CFLAGS) -o hmmsw  hmmsw.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmmt: 	$(OBJ) $(SQUIDOBJ) train_main.o
	$(CC) $(CFLAGS) -o hmmt train_main.o $(OBJ) $(SQUIDOBJ) $(LIBS)

hmm_convert: $(OBJ) $(SQUIDOBJ) convert_main.o
	$(CC) $(CFLAGS) -o hmm_convert convert_main.o $(OBJ) $(SQUIDOBJ) $(LIBS)

align_homologues: $(OBJ) $(SQUIDOBJ) align_homologues.o
	$(CC) $(CFLAGS) -o align_homologues align_homologues.o $(OBJ) $(SQUIDOBJ) $(LIBS)

test:
	(cd Testsuite; perl Run_Tests)

install: 
	test -d $(BINDIR) || mkdir -p $(BINDIR)
	test -d $(MANDIR)/man$(MANSUFFIX) || mkdir -p $(MANDIR)/man$(MANSUFFIX)
	cp $(PROGS) $(BINDIR)/
	for manfile in $(MANS); do\
	  $(INSTMAN) Man/$$manfile.man $(MANDIR)/man$(MANSUFFIX)/$$manfile.$(MANSUFFIX);\
	done

bindist:
	cp $(PROGS) ../hmmer-$(RELEASE)-bin/

clean:
	-rm -f *.o *~ Makefile.bak core $(PROGS) TAGS

tags:
	etags -t $(SRC) $(HDRS)

lint:
	lint $(SRC) $(SQUIDSRC) -lm 

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -c $<		

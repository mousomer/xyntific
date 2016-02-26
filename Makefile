CC 	= gcc
CFLAGS 	= -std=c99 -Wall -lm
OPTIMIZING = -O3 -fPIC -fwrapv
DEBUG = -O0 -ggdb -g3 -Wextra -DDEBUG

LIBS 	= xynMath.h xynDisp.h xynWave.h xynMem.h xynRawVid.h libxyn.h

SRCS 	= xynMem.c xynString.c xynParse.c xynMath.c xynDisp.c xynWave.c genRandFile.c xynRawVid.c\
                displayArray.c displayWavelet.c dwt1d.c idwt1d.c dwt2d.c idwt2d.c calcMoment.c\
                displayVidPosFile.c test.c txt2bin.c

OBJECTS	= $(SRCS:.c=.o)

EXECUTABLES	= genRandFile genLinearFile displayArray displayWavelet \
                 dwt1d idwt1d dwt2d idwt2d dwt2d_buffs dwt2d_buffs_d calcMoment test1 displayVidPosFile txt2bin
EXECOBJ = $(EXECUTABLES:=.o)


all: CFLAGS += $(OPTIMIZING)
all:	$(EXECUTABLES)

debug: CFLAGS += $(DEBUG)
debug:	$(EXECUTABLES)

displayVidPosFile: displayVidPosFile.o xynParse.o xynString.o xynRawVid.o xynDisp.o xynMath.o
	$(CC) $(CFLAGS) -o $@ $^

test1: test1.o xynParse.o xynDisp.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

txt2bin: txt2bin.o xynParse.o xynDisp.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

genLinearFile: genLinearFile.o xynParse.o xynDisp.o xynString.o xynMath.o
	$(CC) $(CFLAGS) -o $@ $^

displayWavelet: displayWavelet.o xynParse.o xynString.o xynWave.o xynMath.o xynDisp.o
	$(CC) $(CFLAGS) -o $@ $^

displayArray: displayArray.o xynParse.o xynString.o xynDisp.o
	$(CC) $(CFLAGS) -o $@ $^

genRandFile: genRandFile.o xynParse.o xynDisp.o xynString.o xynMath.o
	$(CC) $(CFLAGS) -o $@ $^

dwt1d: dwt1d.o  xynParse.o xynWave.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

dwt2d: dwt2d.o  xynParse.o xynWave.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

idwt1d: idwt1d.o xynParse.o xynWave.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

idwt2d: idwt2d.o  xynParse.o xynWave.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

dwt2d_buffs: dwt2d_buffs.o  xynParse.o xynWave.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

dwt2d_buffs_d: dwt2d_buffs_d.o  xynParse.o xynWave.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^

calcMoment: calcMoment.o xynParse.o xynDisp.o xynMath.o xynString.o
	$(CC) $(CFLAGS) -o $@ $^


.o: $(@:.o=.c) libxyn.h
	echo file $< dependencies $^
	$(CC) $(CFLAGS) -c $<

.d: .c
	$(CC) -MM $< -o $@


clean: 
	rm -f $(OBJECTS) $(EXECUTABLES) $(EXECOBJ)


# End of makefile
# DO NOT DELETE

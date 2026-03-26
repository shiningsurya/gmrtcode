
CC     := gcc -DTIMING
CFLAGS := -g -Wall -O2 -Warray-bounds=2
LFLAGS := -lm -lcfitsio -lfftw3

.PHONY: clean

TARGETS := singlepol8bit dualpol8bit 

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

singlepol8bit: gmrtfits.o singlepol8bit.o
	$(CC) $(LFLAGS) -o $@ $^

dualpol8bit: gmrtfits.o dualpol8bit.o
	$(CC) $(LFLAGS) -o $@ $^

clean:
	@rm -fr *.o
	@rm -fr $(TARGETS)

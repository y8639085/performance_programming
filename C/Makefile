SRC=MD.c control.c
OBJ=$(SRC:.c=.o)
CC=icc 
CFLAGS= -g -Ofast -restrict -qopenmp-simd -ipo -xCORE-AVX2
REPORTFLAGS= -qopt-report-file=compilerreport.txt -qopt-report=5

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

MD: $(OBJ)
	$(CC) $(CFLAGS) $(REPORTFLAGS) -o $@ $(OBJ) -lm

all:MD

clean:
	rm -f MD $(OBJ)ls control.o
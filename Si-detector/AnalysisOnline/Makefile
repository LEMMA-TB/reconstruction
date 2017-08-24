CC=c++
CFLAGS=-g -I. $(shell root-config --cflags) 
LIBS=-g $(shell root-config --libs) 

all: AnalyseData

AnalyseData.o: AnalyseData.C
	$(CC) $^ -o $@ -c $(CFLAGS)
AnalyseData: AnalyseData.o
	$(CC) $^ -o $@ $(LIBS)

clean:
	rm AnalyseData
	rm -f *.o

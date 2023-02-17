# Useful flags
DEBUG = -g -Wall -fsanitize=address -fsanitize=leak -fsanitize=undefined
VAL = --track-origins=yes --leak-check=full
CACH = --tool=cachegrind
C_STD = -std=c++14

# Dependecies
DEPS = deps/branched_flow.cpp
DEPS_O = deps/inter/interpolate.o deps/random64.o deps/file_handler.o
INTER = deps/inter/optimization.cpp deps/inter/integration.cpp deps/inter/statistics.cpp deps/inter/ap.cpp deps/inter/alglibmisc.cpp deps/inter/interpolation.cpp deps/inter/linalg.cpp deps/inter/dataanalysis.cpp deps/inter/alglibinternal.cpp deps/inter/specialfunctions.cpp deps/inter/solvers.cpp
INTER_O =  optimization.o integration.o statistics.o ap.o alglibmisc.o interpolation.o linalg.o dataanalysis.o alglibinternal.o specialfunctions.o solvers.o

all : main.o deps/inter/interpolate.o deps/random64.o deps/file_handler.o

main.o : deps/inter/interpolate.o deps/random64.o deps/file_handler.o main.cpp $(DEPS)
	g++ $(C_STD) -O2 deps/inter/interpolate.o deps/random64.o deps/file_handler.o main.cpp $(DEPS) -o $@
	time ./main.o

deps/inter/interpolate.o :
	time g++ $(C_STD) -O2 -c $(INTER)
	ld -r $(INTER_O) -o deps/inter/interpolate.o
	rm -f $(INTER_O)

deps/random64.o : 
	g++ $(C_STD) -O2 -c deps/random64.cpp
	mv random64.o deps/random64.o

deps/file_handler.o :
	g++ $(C_STD) -O2 -c deps/file_handler.cpp
	mv file_handler.o deps/file_handler.o

.PHONY : assembler
assembler : main.cpp $(DEPS) *.h
	gcc $(C_STD) -g -S $(DEPS_O) main.cpp $(DEPS) -o main.S

.PHONY : parallel
parallel : main.cpp
	g++ $(C_STD) -fopenmp -O2 $< $(DEPS) -o parallel.o
	time ./parallel.o

.PHONY : debug
debug : main.cpp
	g++ $(C_STD) $(DEBUG) $(DEPS_O) $< $(DEPS) -o debug.o
	./debug.o

.PHONY : valgrind
valgrind : main.cpp
	g++ $(C_STD) -g $(DEPS_O) $< $(DEPS) -o valgrind.o
	valgrind $(VAL) ./valgrind.o

.PHONY : cachegrind
cachegrind : main.cpp
	g++ $(C_STD) -g $(DEPS_O) $< $(DEPS) -o cachegrind.o
	valgrind $(CACH) ./cachegrind.o

.PHONY : gprof
gprof : main.cpp
	g++ $(C_STD) -Wall -pg $(DEPS_O) $< $(DEPS) -o gprof.o
	./gprof.o
	gprof gprof.o gmon.out > analysis.txt
	most analysis.txt

.PHONY : perf
perf : main.cpp
	g++ $(C_STD) -Wall -pg $(DEPS_O) $< $(DEPS) -o perf.o
	perf record ./perf.o ; perf report

.PHONY : clean
clean :
	rm -f *.o *.txt *.out* debug main *.data



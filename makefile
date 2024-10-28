RM = rm -rf
O_SRCS = tsp.o 
CPP_SRCS = MyThread.cpp main.cpp threads.cpp tsp.cpp twoOpt.cpp direct_flat_graph.cpp dsu.cpp par_boruvka_openmp.cpp
OBJS = MyThread.o main.o threads.o tsp.o twoOpt.o direct_flat_graph.o dsu.o par_boruvka_openmp.o
CPP_DEPS = MyThread.d main.d threads.d tsp.d twoOpt.d direct_flat_graph.d dsu.d par_boruvka_openmp.d
LIBS = -lm -L/usr/include -lpthread -fopenmp


%.o: %.cpp
	mpicxx -g -Wall -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"

tsp: $(OBJS)
	mpicxx  -o "tsp" $(OBJS) $(LIBS) 

all: tsp





clean:
	rm tsp *.o *.d
	rm -f *.tour


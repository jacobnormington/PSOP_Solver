CC = g++
CG = gcc
LINK = -o
CXX_VERSION = -std=c++14
COMPILE = -c
OPTIMIZATION = -O3
CXXFLAG = -Wall $(CXX_VERSION)
SRC = ./src
PBB_LIB = ./lib

OBJDIR = obj

PBB_OBJS = $(patsubst %,$(OBJDIR)/%,$(_PBB_OBJS))

_PBB_OBJS = main.o solver.o hash.o hungarian.o history.o

PROG = sop_solver

$(PROG): $(PBB_OBJS)
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $(LINK) $(PROG) $^

$(OBJDIR)/%.o: $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $(COMPILE) $< -o $@ 

$(OBJDIR)/%.o: $(PBB_LIB)/%.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $(COMPILE) $< -o $@

$(OBJDIR)/%.o: $(LKH_LIB)/%.c
	$(CG) $(OPTIMIZATION) $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o sop_solver
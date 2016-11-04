# Various flags
#CXX  = clang++
CXX = g++-5
LINK = $(CXX)
#CXXFLAGS = -std=c++11 -Wall -g 
CXXFLAGS = -std=c++11 -Wall -O3
#CXXFLAGS = -std=c++11 -Wall -O3 -fopenmp
LFLAGS = -lm
#LFLAGS = -lm -fopenmp
# ifneq "$(findstring noomp, $(MAKECMDGOALS))" "noomp"
# 	CXXFLAGS += -fopenmp
# 	LFLAGS += -fopenmp
# endif

TARGET  = wjs-kmedoids++-lumping

HEADER  = wjs-kmedoids++-lumping.h
FILES = wjs-kmedoids++-lumping.cc

OBJECTS = $(FILES:.cc=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $^ $(LFLAGS) -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




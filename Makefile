# Various flags
CXX  = clang++
LINK = $(CXX)
# CXXFLAGS = -std=c++11 -Wall -g
CXXFLAGS = -std=c++11 -Wall -O3
LFLAGS = -lm

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





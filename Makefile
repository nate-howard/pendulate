#
# compiler
#
CC        = g++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
CXXFLAGS 		= -Wall -Wextra -std=c++11 

#
# scots 
#
SCOTSROOT		= /home/nate/scots
SCOTSINC		= -I$(SCOTSROOT)/utils -I$(SCOTSROOT)/bdd 

#
# cudd 
#
CUDDPATH		=  $(SCOTSROOT)/cudd-3.0.0
CUDDINC 		= -I$(CUDDPATH)
CUDDLIBS		= -lcudd 
CUDDLPATH   	= -L$(CUDDPATH)/lib

BUILDDIR		= ./build

TARGET = pendulate

all: $(TARGET)

%.o:%.cc
	$(CC) -c $(CXXFLAGS) $(CUDDINC) $(SCOTSINC) $< -o $(BUILDDIR)/$@

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) -o $(BUILDDIR)/$(TARGET) $(BUILDDIR)/$(TARGET).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm  ./$(TARGET)  ./$(TARGET).o

run:
	./build/$(TARGET)

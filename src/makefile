
.SUFFIXS 	: .cc .o
CC              = g++ 
LIB             = -L../lib/ -lAlgLib -lbiofiles -lutil -pthread -lstatisticaltest -L/usr/local/lib/ -lz 
CFLAGS          = -O2 -g -Wall -std=c++0x 

SRCS		= main.cc	preprocess.cc	util.cc	calcInfo.cc	#report.cc	
OBJECTS     = main.o	preprocess.o	util.o	calcInfo.o	#report.o
TARGET 		= varinfo 



all: $(TARGET)

$(TARGET) : $(OBJECTS)
		$(CC) -o ../$(TARGET) $(CFLAGS) $(OBJECTS) $(LIB)


dep :
		$(CC) -M $(SRCS) > .depend

clean:
	rm -f $(OBJECTS) .depend
cleanall: clean
	rm -f $(OBJECTS) ../bin/$(TARGET) .depend


ifeq (.depend,$(wildcard .depend))
include .depend
endif




# NOTE: By default, when you type 'make' command and do not provide a file,
#	then the command will automatically use the 'Makefile' as the file. If you
#	have multiple 'Makefile's in a directory, then you have to do:
#		make -f 'name_of_file'

# The previous make file is inefficient for large projects, so
#	we will want to utilize variables, patterns, and dependency
# tracking in order to create a make file that will be more
# easily extendable for larger projects

# Define variables
#		':=' evaluates only once, so it is what it is now
#		'=' evaluates when needed, so if something it depends on changes, it also changes
CC := gcc
CXX := g++
CFLAGS := -Wall -Werror
CXXFLAGS := -Wall -Werror -std=c++11

# This is just telling what gets made when calling 'make' or 'make all'
all: output

# To use variables, you encapsulate the variable name inside $(...)
output: main.o hello.o
	$(CXX) $(CXXFLAGS) main.o hello.o -o output.app

main.o: main.cpp hello.h
	$(CXX) $(CXXFLAGS) -c main.cpp

hello.o: hello.cpp hello.h
	$(CXX) $(CXXFLAGS) -c hello.cpp

clean:
	rm *.o output.app

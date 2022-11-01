# NOTE: By default, when you type 'make' command and do not provide a file,
#	then the command will automatically use the 'Makefile' as the file. If you
#	have multiple 'Makefile's in a directory, then you have to do:
#		make -f 'name_of_file'



# NOTE: Need to figure out how to compile .c and .cpp files
#		(just use same compiler)

# Now we will create a more organized project folder where
#	cpp files go into a 'src' folder, headers into a 'inc' folder,
#	and objects/deps go into a 'bin' folder

# Define names of the folder variables
# NOTE: Add list of all folders here
SRC := src src/adder src/vec
INC := inc inc/adder inc/vec
BIN := bin

# Defined message colors
RED = \033[1;31m
GREEN = \033[1;32m
BLUE = \033[1;34m
YELLOW = \033[1;33m
NC = \033[1;0m

# Cannot color text in windows
ifeq      ($(TARGET_OS),mingw32)
else ifeq ($(TARGET_OS),cygwin)
else
RED = ""
GREEN = ""
BLUE = ""
YELLOW = ""
NC = ""
endif

# Define variables
#		':=' evaluates only once, so it is what it is now
#		'=' evaluates when needed, so if something it depends on changes, it also changes
CC := gcc
CXX := g++
CFLAGS := -O3 -Wall -Werror
CXXFLAGS := -O3 -Wall -Werror -std=c++14

# Additional variables to make the target name and list of objects
# 	Use 'wildcard' function to go into the 'src' directory and create
#		a list of all .cpp files with the directory name $(SRC) appended
#		to the front.
#
#	Then, change 'objs' to have the same format as 'source' (i.e., without the
#		dir), but replace all .cpp with .o. To remove the dir, we use the
#		'notdir' function.
#
#	Then, that is also embedded into the patsubst function to add the bin
#		directory to the objects (which is where we want them)
target := run

c_source := $(wildcard *.c $(foreach fd, $(SRC), $(fd)/*.c))
cpp_source := $(wildcard *.cpp $(foreach fd, $(SRC), $(fd)/*.cpp))
source := $(cpp_source) $(c_source)

# Search top level down for all items (can use wildcard)
SRC_DIRS ?= ./
ex := $(shell find $(SRC) -name 'adder_wrapper.cpp')

c_objs := $(patsubst %, $(BIN)/%, $(notdir $(c_source:.c=.o) ) )
cpp_objs := $(patsubst %, $(BIN)/%, $(notdir $(cpp_source:.cpp=.o) ) )
objs := $(cpp_objs) $(c_objs)

# Create list of include directories from INC var
idirs := $(addprefix -I, $(INC))

# This is just telling what gets made when calling 'make' or 'make all'
# 	In this case, we specify that the 'dir' for 'bin' must be created
#		prior to building our target
all: dir $(target)

# The 'dir' recipie just issues a command line instruction to make a
#		directory named 'bin'. The -p argument prevents errors from occuring
#		if such a directory already exists.
dir:
	mkdir -p bin

# To use variables, you encapsulate the variable name inside $(...)
#		Automatic variables:
#				$^ = replaces name of all dependencies
#				$@ = replaces the name of the target
#				$< = replaces name of the first dependency
run: $(objs)
	@echo "$(RED)Linking objects...$(NC)"
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "$(BLUE)Finished!$(NC)"

# New variables for dependency tracking
#   patsubst = Function called pattern substitution
#         SYNTAX:   $(patsubst 'pattern', 'replacement', 'text')
deps := $(patsubst %.o,%.d,$(objs))
-include $(deps)
DEPFLAGS = -MMD -MF $(@:.o=.d)

# Use % as a pattern similar to *
# NOTE: We have now added a -I option here for the includes
#			This also required us to specify the actual output with
#			-o $@, otherwise, it would put the objects in the wrong folder

# NOTE: You need to specify each 'dependency' for the 'target'
#			with its associated folder. So for every folder in src,
#		  there are a set of instructions
$(BIN)/%.o: src/adder/%.cpp
	@echo "$(GREEN)Compiling...$(NC)"
	$(CXX) $(CXXFLAGS) -c $< $(DEPFLAGS) $(idirs) -o $@

$(BIN)/%.o: src/vec/%.cpp
	@echo "$(GREEN)Compiling...$(NC)"
	$(CXX) $(CXXFLAGS) -c $< $(DEPFLAGS) $(idirs) -o $@

$(BIN)/%.o: src/%.cpp
	@echo "$(GREEN)Compiling...$(NC)"
	$(CXX) $(CXXFLAGS) -c $< $(DEPFLAGS) $(idirs) -o $@

# How to handle .c files
#$(BIN)/%.o: $(SRC)/%.c
#	@echo "$(GREEN)Compiling...$(NC)"
#	$(CXX) $(CFLAGS) -c $< $(DEPFLAGS) $(idirs) -o $@

# Use the .PHONY target as a list of names to append to the make command
.PHONY: clean help

clean:
	rm $(objs) $(deps) $(target)

help:
	@echo "SRC: $(SRC)"
	@echo "cpp_source: $(cpp_source)"
	@echo "cpp_objs: $(cpp_objs)"
	@echo "source: $(source)"
	@echo "objs: $(objs)"
	@echo "deps: $(deps)"
	@echo "idirs: $(idirs)"

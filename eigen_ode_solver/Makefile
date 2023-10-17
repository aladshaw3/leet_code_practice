# NOTE: By default, when you type 'make' command and do not provide a file,
#	then the command will automatically use the 'Makefile' as the file. If you
#	have multiple 'Makefile's in a directory, then you have to do:
#		make -f 'name_of_file'

# Define names of the folder variables
# NOTE: Add list of all folders here
SRC := src 
INC := inc ../../eigen
BIN := bin

# Define standard console colors
BLACK        := $(shell tput -Txterm setaf 0)
RED          := $(shell tput -Txterm setaf 1)
GREEN        := $(shell tput -Txterm setaf 2)
YELLOW       := $(shell tput -Txterm setaf 3)
LIGHTPURPLE  := $(shell tput -Txterm setaf 4)
PURPLE       := $(shell tput -Txterm setaf 5)
BLUE         := $(shell tput -Txterm setaf 6)
WHITE        := $(shell tput -Txterm setaf 7)
NC := $(shell tput -Txterm sgr0)

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
target := run_sim

c_source := $(wildcard *.c $(foreach fd, $(SRC), $(fd)/*.c))
cpp_source := $(wildcard *.cpp $(foreach fd, $(SRC), $(fd)/*.cpp))
source := $(cpp_source) $(c_source)

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
	@echo "$(BLUE)Creating bin/...$(NC)"
	mkdir -p bin

# To use variables, you encapsulate the variable name inside $(...)
#		Automatic variables:
#				$^ = replaces name of all dependencies
#				$@ = replaces the name of the target
#				$< = replaces name of the first dependency
$(target): $(objs)
	@echo "$(YELLOW)Linking objects...$(NC)"
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "$(GREEN)Finished!$(NC)"

# Variables for dependency tracking
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

$(BIN)/%.o: src/%.cpp
	@echo "$(RED)Compiling...$(NC)"
	$(CXX) $(CXXFLAGS) -c $< $(DEPFLAGS) $(idirs) -o $@

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

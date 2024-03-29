# NOTE: Need to figure out how to compile .c and .cpp files
#		(just use same compiler)

# Now we will create a more organized project folder where
#	cpp files go into a 'src' folder, headers into a 'inc' folder,
#	and objects/deps go into a 'bin' folder

# Define names of the folder variables
SRC := src
INC := inc
BIN := bin

# Defined message colors
RED = \033[1;31m
GREEN = \033[1;32m
BLUE = \033[1;34m
YELLOW = \033[1;33m
NC = \033[1;0m

# Define variables
#		':=' evaluates only once, so it is what it is now
#		'=' evaluates when needed, so if something it depends on changes, it also changes
CC := gcc
CXX := g++
CFLAGS := -O3 -Wall -Werror
CXXFLAGS := -O3 -Wall -Werror -std=c++11

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
target := output

c_source := $(wildcard $(SRC)/*.c)
cpp_source := $(wildcard $(SRC)/*.cpp)
source := $(cpp_source) $(c_source)

c_objs := $(patsubst %, $(BIN)/%, $(notdir $(c_source:.c=.o) ) )
cpp_objs := $(patsubst %, $(BIN)/%, $(notdir $(cpp_source:.cpp=.o) ) )
objs := $(cpp_objs) $(c_objs)

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
output: $(objs)
	@echo "$(RED)Linking objects...$(NC)"
	$(CXX) $(CXXFLAGS) $^ -o $@.app
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
$(BIN)/%.o: $(SRC)/%.cpp
	@echo "$(GREEN)Compiling...$(NC)"
	$(CXX) $(CXXFLAGS) -c $< $(DEPFLAGS) -I$(INC) -o $@

# How to handle .c files
$(BIN)/%.o: $(SRC)/%.c
	@echo "$(GREEN)Compiling...$(NC)"
	$(CXX) $(CFLAGS) -c $< $(DEPFLAGS) -I$(INC) -o $@

# Use the .PHONY target as a list of names to append to the make command
.PHONY: clean help

clean:
	rm $(objs) $(deps) $(target).app

help:
	@echo "source: $(source)"
	@echo "objs: $(objs)"
	@echo "deps: $(deps)"

# NOTE: We are still assuming that all files are in same directory

# In this makefile, we will add console messages and colors to the commands
#	to better track the progress of the build.
#
#	Adding messages is accomplished through the @echo command
# Adding colors is done through special color info
#		NOTE: NC = no color
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
CFLAGS := -Wall -Werror
CXXFLAGS := -Wall -Werror -std=c++11

# Additional variables to make the target name and list of objects
target := output
objs := main.o hello.o

# This is just telling what gets made when calling 'make' or 'make all'
all: $(target)

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
%.o: %.cpp
	@echo "$(GREEN)Compiling...$(NC)"
	$(CXX) $(CXXFLAGS) -c $< $(DEPFLAGS)

# Use the .PHONY target as a list of names to append to the make command
.PHONY: clean help

clean:
	rm $(objs) $(deps) $(target).app

help:
	@echo "objs: $(objs)"
	@echo "deps: $(deps)"

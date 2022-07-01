# To compile main manually: g++ main.cpp hello.cpp

# Instead, we can create a formula to automate this

# That formula follows this basic recipie:
#Target: dependencies
#		command
#		(another command)...

# NOTE: Running 'make' with no additial argument will always run
#	the FIRST commant set provided (in this case, 'output')


# So, we are creating an executable target named 'output'
#	which depends on the objects of 'main.o' and 'hello.o'
#	and is created with the command: g++ main.o hello.o -o output
#	This command is using g++ to link objects main and hello together
#	to create an executable named output
output: main.o hello.o
	g++ main.o hello.o -o output

# Now we should define how main.o and hello.o are created
#	The '-c' flag instructs the build to just compile the
#	main.cpp into and object, without building an executable

# NOTE: This instruction is only run if there is a change in
#	main.cpp. However, this may not be the behavior we want, because
#	main.cpp includes the header 'hello.h'. Thus, if we change something
# in that header, this may not get rebuilt. 
main.o: main.cpp
	g++ -c main.cpp

hello.o: hello.cpp hello.h
	g++ -c hello.cpp

# This is a target witout dependencies
#	it will run the commands given below to remove
#	the object files and executable
clean:
	rm *.o
	rm output

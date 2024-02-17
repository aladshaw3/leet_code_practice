// Include the hello functions
#include "hello.h"

// The main function is the primary executable
int main()
{
    hello();
	Algo<std::chrono::milliseconds, 500> test;
	auto in = std::make_shared<float>(1.0);
	std::cout << in.use_count() << std::endl;
	auto out = std::make_shared<float>(0.0);
	
	std::cout << *in << std::endl;
	*in = 5;
	std::cout << *in << std::endl;
	
	std::cout << *out << std::endl;
	test.Compute(*in, *out);
	std::cout << *out << std::endl;
	
	//std::cout << in.use_count() << std::endl;
	
	test.Compute(in, out);
	std::cout << *out << std::endl;
	
	test.Compute(in, out);
	std::cout << *out << std::endl;
    return 0;
}

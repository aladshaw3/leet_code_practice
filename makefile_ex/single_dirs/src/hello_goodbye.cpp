// Add include for the "hello_goodbye.h"
#include "hello_goodbye.h"

// Define how the function 'hello' behaves
void hello_goodbye()
{
  // Call the function from hello.h, which was included in hello_goodbye.h
  hello();

  // Add an additional message
  std::cout << "Goodbye" << std::endl;
}

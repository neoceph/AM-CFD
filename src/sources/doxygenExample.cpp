#include "doxygenExample.hpp"


/// Say hello
/// @param message The message to print
/* this is an example text for testing Doxygen
    The idea is to check whether documentation works
*/
void Foo::say_hello(std::string message) const {
  std::cout << "Hello: " << message << std::endl;
}

#ifndef DOXYGENEXAMPLE_HPP
#define DOXYGENEXAMPLE_HPP

#include <iostream>
#include <string>
class Foo {

  /// Say hello
  /// @param message The message to print
  /* this is an example text for testing Doxygen
      The idea is to check whether documentation works
  */
  void say_hello(std::string message) const;


  void say_helloAgain(std::string message, std::string name) const;
};

#endif /* DOXYGENEXAMPLE_HPP */
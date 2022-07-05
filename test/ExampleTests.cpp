#include "Example.hpp"
#include<gtest/gtest.h>

TEST(ExampleTests, MAC){

    int x = 42;
    int y = 16;
    int sum = 100;
    int oldSum = sum;
    int expectedNewSum = oldSum + x * y;

    EXPECT_EQ(
        expectedNewSum,
        MAC(x, y, sum)
        );
    EXPECT_EQ(
        expectedNewSum,
        sum
    )<< "The sum is not equal to the expected new sum!";
    

}

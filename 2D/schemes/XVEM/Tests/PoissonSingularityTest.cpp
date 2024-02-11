#include "PoissonSingularityTest.h"

namespace PolyMesh2D
{
    TEST_F(PoissonSingularityTest, TestWithDefaultInput)
    {
        // Your test logic here
        // Example: Testing PoissonSingularity function
        const Functional::Function<2, 1> result(PoissonSingularity(1.5));

        ASSERT_EQ(4, 4);
        
        // Add assertions to validate the output based on expected behavior
        ASSERT_NEAR(result.value(Eigen::Vector2d(0.5, 0.3)), 0.32268588784198498, 1E-15);
        ASSERT_NEAR(result.derivative(Eigen::Vector2d(0.5, 0.3))(0), 0.30574840665177031, 1E-15);
        ASSERT_NEAR(result.derivative(Eigen::Vector2d(0.5, 0.3))(1), 1.1038487614569744, 1E-15);
        // Add more assertions as needed
    }
}

// Define the main function to run the tests
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
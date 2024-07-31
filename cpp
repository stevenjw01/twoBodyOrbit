#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>
#include "OSQP.h"
#include "nml.h"
#include "OSQP_print.h"
#include "Gruven.h"
#include "nml_util.h"

// Define extern C functions
extern "C" {
    struct osqp_in osqp_nasa_unit_test_main(const char* test);
    struct dimensions get_dimensions(const char *test);
}

// Function to allocate 2D array dynamically
double** alloc_2d(int nrows, int ncols) {
    double** array2d = (double**)calloc(nrows, sizeof(double*));
    for (int i = 0; i < nrows; i++) {
        array2d[i] = (double*)calloc(ncols, sizeof(double));
    }
    return array2d;
}

// Test Fixture for parameterized tests
class QPTest : public ::testing::TestWithParam<std::string> {};

TEST_P(QPTest, OutputWithinRange) {
    std::string testName = GetParam();
    struct osqp_in osqp_in_vect;
    struct osqp_out osqp_out_vect;
    struct dimensions dims;

    osqp_in_vect = osqp_nasa_unit_test_main(testName.c_str());
    dims = get_dimensions(testName.c_str());

    // Adjust as necessary to handle different accuracy levels if required
    osqp_in_vect.rprim_eps = 1e-8;

    // Perform the test
    osqp_out_vect = osqp_lite(osqp_in_vect);

    // Validate output (example: assuming the output is validated against a known good value)
    double expected_output = /* define the expected value based on test name */;
    double tolerance = 0.1;

    double **test_dx = alloc_2d(dims.N, dims.M);
    char newFile[150];
    snprintf(newFile, sizeof(newFile), "%s/%s-Dx-%.3e.csv", FILEPATH, testName.c_str(), osqp_in_vect.rprim_eps);

    FILE *qp_test = fopen(newFile, "r");
    if (qp_test == NULL) {
        fprintf(stderr, "\n Error reading test_data reference file! \n");
        return;
    }

    for (int ii = 0; ii < dims.M; ii++) {
        for (int jj = 0; jj < dims.N; jj++) {
            fscanf(qp_test, "%lf,", &test_dx[ii][jj]);
        }
    }

    fclose(qp_test);

    for (int ii = 0; ii < dims.M; ii++) {
        for (int jj = 0; jj < dims.N; jj++) {
            EXPECT_NEAR(test_dx[ii][jj], osqp_out_vect.D_x->data[ii][jj], tolerance);
        }
    }
}

// Function to get test names dynamically
std::vector<std::string> getTestNames() {
    // Add the list of your test cases here
    return {"Govind", "Test2", "Test3"};
}

// Instantiate the test case with different parameters
INSTANTIATE_TEST_SUITE_P(DynamicTests, QPTest, ::testing::ValuesIn(getTestNames()));

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

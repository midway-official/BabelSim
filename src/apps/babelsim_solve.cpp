#include "babelsim/application.h"
#include <iostream>

int main(int argc, char* argv[]) {
    return babelsim::runApplication(argc, argv, [](const char* message) {
        std::cerr << "babelsim-solve: " << message << '\n';
    });
}

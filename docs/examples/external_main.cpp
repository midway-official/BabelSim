#include "babelsim/application.h"
#include <iostream>

namespace {
void reportError(const char* message) {
    std::cerr << "tutorial-solver: " << message << '\n';
}
}

int main(int argc, char** argv) {
    return babelsim::runApplication(argc, argv, reportError);
}

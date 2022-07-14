#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
#include <stdio.h>
#include <chrono>

using namespace std;

int main(int argc, char*argv[]) {
    solver s;
    if (argc > 3) {
        cout << "Usage: <instance name> <time limit>";
        exit(-1);
    }
    s.solve(argv[1],atoi(argv[2]));
    return 0;
}
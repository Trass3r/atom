// Copyright 2015 Adobe Systems Incorporated
// All Rights Reserved.

#include "Halide.h"

Halide::Func getFunction();

int main(int argc, char **argv) {
    Halide::Func theFunc = getFunction();

    if (argc >= 3) {
        std::vector<Halide::Argument> arguments = theFunc.infer_arguments();

        Halide::Target target = Halide::get_target_from_environment();
        target.set_feature(Halide::Target::Feature::UserContext);

        theFunc.compile_to_file(argv[1], arguments, argv[2], target);

        return 0;
    }

    return 1;
}

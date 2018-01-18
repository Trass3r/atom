// Copyright 2015 Adobe Systems Incorporated
// All Rights Reserved.

#include <Halide.h>
#include <string>

using namespace Halide;

struct MyFirstGenerator : public Halide::Generator<MyFirstGenerator>
{
	Var x{"x"}, y{"y"}, c{"c"};

	Input<int> size{"size", 32, 4, 256};
	Input<int> shade_1{"shade1", 192, 0, 255};
	Input<int> shade_2{"shade2", 64, 0, 255};

	Output<Buffer<uint8_t>> result{"result", 3};

	void generate()
	{
		result(x, y, c) = cast<uint8_t>(select((floor (x / size) % 2) != (floor(y / size) % 2), shade_1, shade_2));
	}

	void schedule()
	{
		if (auto_schedule)
		{
			result.estimate(c, 0, 3)
				.estimate(x, 0, 512)
				.estimate(y, 0, 512);
		}
		else
		{
			result.parallel(y).vectorize(x, 8);
		}
	}
};

HALIDE_REGISTER_GENERATOR(MyFirstGenerator, generator)

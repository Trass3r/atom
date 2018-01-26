#include <Halide.h>
using namespace Halide;

#define _(x) x{#x}

using namespace Halide::ConciseCasts;

// Shared variables
Var x("x"), y("y"), c("c"), yi("yi"), yo("yo"), yii("yii"), xi("xi");

Var x_i("x_i");
Var x_i_vi("x_i_vi");
Var x_i_vo("x_i_vo");
Var x_o("x_o");
Var x_vi("x_vi");
Var x_vo("x_vo");

// Average two positive values rounding up
Expr avg(Expr a, Expr b)
{
    Type wider = a.type().with_bits(a.type().bits() * 2);
    return cast(a.type(), (cast(wider, a) + b + 1)/2);
}

Expr blur121(Expr a, Expr b, Expr c)
{
    return avg(avg(a, c), b);
}

Func interleave_x(Func a, Func b)
{
    Func out("interleave_x");
    out(x, y) = select((x%2)==0, a(x/2, y), b(x/2, y));
    return out;
}

Func interleave_y(Func a, Func b)
{
    Func out("interleave_y");
    out(x, y) = select((y%2)==0, a(x, y/2), b(x, y/2));
    return out;
}

struct Demosaic final : public Halide::Generator<Demosaic>
{
    Input<Func> input{"input", Int(16), 3}; // in RGGB bands
    Output<Func> output{"output", Int(16), 3};

    Input<bool> hmm{"hmm", false};

    Func _(R), _(Gr), _(Gb), _(B);
    Func _(outGb), _(outGr);

    void generate()
    {
        R (x, y) = input(x + 1, y + 1, 0);
        Gb(x, y) = input(x + 1, y + 1, 1);
        Gr(x, y) = input(x + 1, y + 1, 2);
        B (x, y) = input(x + 1, y + 1, 3);

        // TODO: BoundaryConditions::repeat_edge

        // pattern:
        // G B
        // R G

        Func _(Ghb), _(Ghr), _(Gvb), _(Gvr);
        // even rows = blue
        Ghb(x, y) = (Gb(x, y) + Gb(x+1, y)) / 2  + (2 * B(x, y) - B(x-1, y) - B(x+1, y)) / 4;
        // odd rows = red
        Ghr(x, y) = (Gr(x-1, y) + Gr(x, y)) / 2  + (2 * R(x, y) - R(x-1, y) - R(x+1, y)) / 4;

        // even cols = red
        Gvr(x, y) = (Gb(x, y) + Gb(x, y+1)) / 2  + (2 * R(x, y) - R(x, y-1) - R(x, y+1)) / 4;
        // odd cols = blue
        Gvb(x, y) = (Gr(x, y-1) + Gr(x, y)) / 2  + (2 * B(x, y) - B(x, y-1) - B(x, y+1)) / 4;

		// chrominances of reconstructed images
		Func _(Chb), _(Chr), _(Cvb), _(Cvr);
		Chb(x, y) = Ghb(x, y) - B(x, y);
		Chr(x, y) = Ghr(x, y) - R(x, y);
		Cvb(x, y) = Gvb(x, y) - B(x, y);
		Cvr(x, y) = Gvr(x, y) - R(x, y);

        //output(x, y, c) = select(hmm, Chb(x/2, y/2), Cvb(x/2, y/2));
        //return;

		// gradients of chrominances
		Func _(Dhb), _(Dhr), _(Dvb), _(Dvr);
		Dhb(x, y) = abs(Chb(x, y) - Chb(x+1, y));
		Dhr(x, y) = abs(Chr(x, y) - Chr(x+1, y));
		Dvb(x, y) = abs(Cvb(x, y) - Cvb(x, y+1));
		Dvr(x, y) = abs(Cvr(x, y) - Cvr(x, y+1));

		const int C = 3;
		Func _(deltaHb), _(deltaVb), _(deltaHr), _(deltaVr);
		deltaHb(x, y) = C * Dhb(x, y) + C * Dhb(x-1, y) + Dhb(x-1, y-1) + Dhb(x-1, y+1) + Dhb(x, y-1) + Dhb(x, y+1) + Dhr(x, y) + Dhr(x, y+1);
		deltaHr(x, y) = C * Dhr(x, y) + C * Dhr(x-1, y) + Dhr(x-1, y-1) + Dhr(x-1, y+1) + Dhr(x, y-1) + Dhr(x, y+1) + Dhb(x-1, y-1) + Dhb(x-1, y);
		deltaVb(x, y) = C * Dvr(x, y) + C * Dvr(x-1, y) + Dvr(x-1, y-1) + Dvr(x-1, y+1) + Dvr(x, y-1) + Dvr(x, y+1) + Dvb(x-1, y-1) + Dvb(x-1, y);
		deltaVr(x, y) = C * Dvb(x, y) + C * Dvb(x-1, y) + Dvb(x-1, y-1) + Dvb(x-1, y+1) + Dvb(x, y-1) + Dvb(x, y+1) + Dvr(x, y) + Dvr(x, y+1);

        // pattern:
        // G B G B G B G B
        // R G R G R G R G
        // G B G[B]G B G B
        // R G R G R G R G
        // G B G B G B G B
        Expr Eb = deltaVb(x, y) < deltaHb(x, y);
        Expr Er = deltaVr(x, y) < deltaHr(x, y);
        outGb(x, y) = select(Eb, Gvb(x, y), Ghb(x, y));
		outGr(x, y) = select(Er, Gvr(x, y), Ghr(x, y));

        Func _(outG);
        //Func _(Gh), _(Gv);
		// Gh(x, y) = select(y % 2 == 0, Ghb, Ghr);
		// Gv(x, y) = select(x % 2 == 0, Gvr, Gvb);
        outG(x, y) = select(y % 2 == 0,
            select(x % 2 == 0, Gb(x/2, y/2), outGb(x/2, y/2)),
            select(x % 2 == 0, outGr(x/2, y/2), Gr(x/2, y/2))
        );

        Func _(Rgr), _(Rgb), _(Bgr), _(Bgb);
        // G B G B
        // R G R G
        // G[B]G B
        // R G R G

        // green pixels: red through bilinear interpolation of R-G
        Bgb(x, y) = Gb(x, y) + ( B(x-1, y)-Gb(x-1, y) + B(x, y)-Gb(x-1, y) )/2;
        Rgb(x, y) = Gb(x, y) + ( R(x-1, y)-Gr(x-1, y) + R(x, y)-Gr(x-1, y) )/2;
        Bgr(x, y) = Gb(x, y) + ( B(x, y)-Gb(x, y) + B(x+1, y)-Gb(x+1, y) )/2;
        Rgr(x, y) = Gr(x, y) + ( R(x, y)-Gr(x, y) + R(x+1, y)-Gr(x+1, y) )/2;

        // reconstruction of the red and blue values in the blue and red pixels,
        // respectively.

        // R-B = R-(Gb-1) + R-(Gb+1)
        Func _(Rb), _(Br);
        Rb(x, y) = select(Er,
            B(x, y) + (Rgb(x, y) - Bgb(x, y) + Rgb(x + 1, y) - Bgb(x + 1, y)) / 2 // (B(x-1, y) B(x+1, y)) / 2
            ,
            B(x, y) + (Rgb(x, y) - Bgb(x, y) + Rgb(x, y + 1) - Bgb(x, y + 1)) / 2 // (B(x-1, y) B(x+1, y)) / 2
        );
        Br(x, y) = select(Eb,
            R(x, y) + (Bgr(x, y) - Rgr(x, y) + Bgr(x + 1, y) - Rgr(x + 1, y)) / 2 // (B(x-1, y) B(x+1, y)) / 2
            ,
            R(x, y) + (Bgr(x, y) - Rgr(x, y) + Bgr(x, y + 1) - Rgr(x, y + 1)) / 2 // (B(x-1, y) B(x+1, y)) / 2
        );

        //output(x, y, c) = outG(x, y);
        // G B G B
        // R G R G
        output(x, y, c) = select(c == 1,
            outG(x, y),
            c == 0,
            select(y % 2 == 0,
                select(x % 2 == 0, Bgb(x/2, y/2), B(x/2, y/2)),
                select(x % 2 == 0, Br(x/2, y/2), Bgr(x/2, y/2))
            ),
            select(y % 2 == 0,
                select(x % 2 == 0, Rgb(x/2, y/2), Rb(x/2, y/2)),
                select(x % 2 == 0, R(x/2, y/2), Rgr(x/2, y/2))
            )
        );
    }

    void schedule()
    {
        if (auto_schedule)
            return;

        R.compute_root();
        Gr.compute_root();
        Gb.compute_root();
        B.compute_root();
        outGb.compute_root();
        outGr.compute_root();
    }
};

struct CameraPipe final : public Halide::Generator<CameraPipe>
{
    Input<Buffer<uint8_t>> input{"input", 3};
    Output<Buffer<uint8_t>> processed{"processed", 3};

    Input<bool> showRaw{"showRaw", false};
    Input<bool> hmm{"hmm", false};

    void generate();
    void schedule();

    Func _(deinterleaved);
    Func _(base);
    Func _(raw);
    std::unique_ptr<Demosaic> demosaiced;
private:
    Func deinterleave(Func raw);
};

Func CameraPipe::deinterleave(Func raw)
{
    // Deinterleave the color channels
    Func _(deinterleaved);
    deinterleaved(x, y, c) = select(c == 0, raw(2*x, 2*y),
                                    c == 1, raw(2*x+1, 2*y),
                                    c == 2, raw(2*x, 2*y+1),
                                            raw(2*x+1, 2*y+1));
    return deinterleaved;
}

void CameraPipe::generate()
{
    base(x, y, c) = cast<int16_t>(BoundaryConditions::mirror_image(input)(x, y, c));
    //base(x, y, c) = cast<int16_t>(input(x, y, c));

    // turn RGB into bayer AND deinterleave
    deinterleaved(x, y, c) = select(c == 1, base(2*x, 2*y, 1),     // green0
                                    c == 2, base(2*x+1, 2*y+1, 1), // green1
                                    c == 0, base(2*x, 2*y+1, 0),   // red
                                            base(2*x+1, 2*y, 2));  // blue

    raw(x, y) = base(x, y, select((x%2) == (y%2), 1,
                        (x%2)==0 && (y%2) == 1, 0,
                        2));
    demosaiced = create<Demosaic>();
    demosaiced->apply(deinterleaved, hmm);

    Expr forDisplay = select(showRaw, raw(x, y), demosaiced->output(x, y, c));
    processed(x, y, c) = cast<uint8_t>(clamp(forDisplay, 0, 255));
}

void CameraPipe::schedule()
{
    if (auto_schedule)
    {
        input.dim(0).set_bounds_estimate(0, 2592);
        input.dim(1).set_bounds_estimate(0, 1968);
        input.dim(2).set_bounds_estimate(0, 3);

        processed
            .estimate(c, 0, 3)
            .estimate(x, 0, 2592)
            .estimate(y, 0, 1968);
    }
    else
    {
        base.compute_root();

        raw
            .compute_at(processed, x_o)
            .split(x, x_vo, x_vi, 16)
            .vectorize(x_vi);

        Expr out_width = processed.width();
        Expr out_height = processed.height();
        int vec = get_target().natural_vector_size(UInt(16));
        Expr strip_size = 32;
        strip_size = (strip_size / 2) * 2;
        processed
            .compute_root()
            .split(x, x_o, x_i, 64)
            .reorder(x_i, y, c, x_o)
            .split(x_i, x_i_vo, x_i_vi, 32)
            .vectorize(x_i_vi)
            //.parallel(x_o)
        ;

// We can generate slightly better code if we know the splits divide the extent.
        processed
            .bound(c, 0, 3)
//          .bound(x, 0, ((out_width)/(2*vec))*(2*vec))
//          .bound(y, 0, (out_height/strip_size)*strip_size)
        ;
    }
}

HALIDE_REGISTER_GENERATOR(CameraPipe, generator)

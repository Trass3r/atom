#include <Halide.h>
using namespace Halide;

#define _(x) x{#x}

using namespace Halide::ConciseCasts;

// Shared variables
Var x("x"), y("y"), c("c"), yi("yi"), yo("yo"), yii("yii"), xi("xi");

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

class Demosaic : public Halide::Generator<Demosaic>
{
public:
    GeneratorParam<LoopLevel> intermed_compute_at{"intermed_compute_at", LoopLevel::inlined()};
    GeneratorParam<LoopLevel> intermed_store_at{"intermed_store_at", LoopLevel::inlined()};
    GeneratorParam<LoopLevel> output_compute_at{"output_compute_at", LoopLevel::inlined()};

    // Inputs and outputs
    Input<Func> deinterleaved{ "deinterleaved", Int(16), 3 };
    Output<Func> output{ "output", Int(16), 3 };

    // Defines outputs using inputs
    void generate()
    {
        // These are the values we already know from the input
        // x_y = the value of channel x at a site in the input of channel y
        // gb refers to green sites in the blue rows
        // gr refers to green sites in the red rows

        // Give more convenient names to the four channels we know
        Func r_r("r_r"), g_gr("g_gr"), g_gb("g_gb"), b_b("b_b");

        g_gr(x, y) = deinterleaved(x, y, 0);
        r_r(x, y)  = deinterleaved(x, y, 1);
        b_b(x, y)  = deinterleaved(x, y, 2);
        g_gb(x, y) = deinterleaved(x, y, 3);

        // These are the ones we need to interpolate
        Func b_r("b_r"), g_r("g_r"), b_gr("b_gr"), r_gr("r_gr"), b_gb("b_gb"), r_gb("r_gb"), r_b("r_b"), g_b("g_b");

        // First calculate green at the red and blue sites

        // Try interpolating vertically and horizontally. Also compute
        // differences vertically and horizontally. Use interpolation in
        // whichever direction had the smallest difference.
        Expr gv_r  = avg(g_gb(x, y-1), g_gb(x, y));
        Expr gvd_r = absd(g_gb(x, y-1), g_gb(x, y));
        Expr gh_r  = avg(g_gr(x+1, y), g_gr(x, y));
        Expr ghd_r = absd(g_gr(x+1, y), g_gr(x, y));

        g_r(x, y)  = select(ghd_r < gvd_r, gh_r, gv_r);

        Expr gv_b  = avg(g_gr(x, y+1), g_gr(x, y));
        Expr gvd_b = absd(g_gr(x, y+1), g_gr(x, y));
        Expr gh_b  = avg(g_gb(x-1, y), g_gb(x, y));
        Expr ghd_b = absd(g_gb(x-1, y), g_gb(x, y));

        g_b(x, y)  = select(ghd_b < gvd_b, gh_b, gv_b);

        // Next interpolate red at gr by first interpolating, then
        // correcting using the error green would have had if we had
        // interpolated it in the same way (i.e. add the second derivative
        // of the green channel at the same place).
        Expr correction;
        correction = g_gr(x, y) - avg(g_r(x, y), g_r(x-1, y));
        r_gr(x, y) = correction + avg(r_r(x-1, y), r_r(x, y));

        // Do the same for other reds and blues at green sites
        correction = g_gr(x, y) - avg(g_b(x, y), g_b(x, y-1));
        b_gr(x, y) = correction + avg(b_b(x, y), b_b(x, y-1));

        correction = g_gb(x, y) - avg(g_r(x, y), g_r(x, y+1));
        r_gb(x, y) = correction + avg(r_r(x, y), r_r(x, y+1));

        correction = g_gb(x, y) - avg(g_b(x, y), g_b(x+1, y));
        b_gb(x, y) = correction + avg(b_b(x, y), b_b(x+1, y));

        // Now interpolate diagonally to get red at blue and blue at
        // red. Hold onto your hats; this gets really fancy. We do the
        // same thing as for interpolating green where we try both
        // directions (in this case the positive and negative diagonals),
        // and use the one with the lowest absolute difference. But we
        // also use the same trick as interpolating red and blue at green
        // sites - we correct our interpolations using the second
        // derivative of green at the same sites.

        correction = g_b(x, y)  - avg(g_r(x, y), g_r(x-1, y+1));
        Expr rp_b  = correction + avg(r_r(x, y), r_r(x-1, y+1));
        Expr rpd_b = absd(r_r(x, y), r_r(x-1, y+1));

        correction = g_b(x, y)  - avg(g_r(x-1, y), g_r(x, y+1));
        Expr rn_b  = correction + avg(r_r(x-1, y), r_r(x, y+1));
        Expr rnd_b = absd(r_r(x-1, y), r_r(x, y+1));

        r_b(x, y)  = select(rpd_b < rnd_b, rp_b, rn_b);

        // Same thing for blue at red
        correction = g_r(x, y)  - avg(g_b(x, y), g_b(x+1, y-1));
        Expr bp_r  = correction + avg(b_b(x, y), b_b(x+1, y-1));
        Expr bpd_r = absd(b_b(x, y), b_b(x+1, y-1));

        correction = g_r(x, y)  - avg(g_b(x+1, y), g_b(x, y-1));
        Expr bn_r  = correction + avg(b_b(x+1, y), b_b(x, y-1));
        Expr bnd_r = absd(b_b(x+1, y), b_b(x, y-1));

        b_r(x, y)  =  select(bpd_r < bnd_r, bp_r, bn_r);

        // Resulting color channels
        Func _(r), _(g), _(b);

        // Interleave the resulting channels
        r = interleave_y(interleave_x(r_gr, r_r),
                         interleave_x(r_b, r_gb));
        g = interleave_y(interleave_x(g_gr, g_r),
                         interleave_x(g_b, g_gb));
        b = interleave_y(interleave_x(b_gr, b_r),
                         interleave_x(b_b, b_gb));

        output(x, y, c) = select(c == 0, r(x, y),
                                 c == 1, g(x, y),
                                         b(x, y));

        // These are the stencil stages we want to schedule
        // separately. Everything else we'll just inline.
        intermediates.push_back(g_r);
        intermediates.push_back(g_b);
    }

    void schedule() {
        Pipeline p(output);

        if (auto_schedule)
            return;

        int vec = get_target().natural_vector_size(UInt(16));
        for (Func f : intermediates)
        {
            f.compute_at(intermed_compute_at)
                .store_at(intermed_store_at)
                .vectorize(x, 2*vec, TailStrategy::RoundUp)
                .fold_storage(y, 4);
        }
        intermediates[1].compute_with(
            intermediates[0], x,
            {{x, LoopAlignStrategy::AlignStart}, {y, LoopAlignStrategy::AlignStart}});
        output.compute_at(output_compute_at)
            .vectorize(x)
            .unroll(y)
            .reorder(c, x, y)
            .unroll(c);
    }

private:
    // Intermediate stencil stages to schedule
    std::vector<Func> intermediates;
};

struct CameraPipe : public Halide::Generator<CameraPipe>
{
    Input<Buffer<uint8_t>> input{"input", 3};
    Output<Buffer<uint8_t>> processed{"processed", 3};

    Input<bool> showRaw{"showRaw", false};
    Input<uint16_t> zoom{"zoom", 1, 1, 10};

    void generate();
    void schedule();

    Func _(deinterleaved);
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
    Func _(base);
    base(x, y, c) = cast<int16_t>(BoundaryConditions::repeat_edge(input)(x, y, c));

    deinterleaved(x, y, c) = select(c == 0, base(2*x, 2*y, 1),
                                    c == 1, base(2*x+1, 2*y, 2),
                                    c == 2, base(2*x, 2*y+1, 0),
                                            base(2*x+1, 2*y+1, 1));

    Func _(raw);
    raw(x, y) = base(x, y, select((x%2) == (y%2), 1,
                        (x%2)==0 && (y%2) == 1, 0,
                        2));
    //Func deinterleaved = deinterleave(shifted);
    demosaiced = create<Demosaic>();
    demosaiced->apply(deinterleaved);

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
        Expr out_width = processed.width();
        Expr out_height = processed.height();
        // In HVX 128, we need 2 threads to saturate HVX with work,
        //and in HVX 64 we need 4 threads, and on other devices,
        // we might need many threads.
        Expr strip_size;
        strip_size = 32;
        strip_size = (strip_size / 2) * 2;

        int vec = get_target().natural_vector_size(UInt(16));
        processed.compute_root()
        .reorder(c, x, y)
            .split(y, yi, yii, 2, TailStrategy::RoundUp)
            .split(yi, yo, yi, strip_size / 2)
            .vectorize(x, 2*vec, TailStrategy::RoundUp)
            .unroll(c)
            .parallel(yo);

        deinterleaved.compute_at(processed, yi).store_at(processed, yo)
            .fold_storage(y, 8)
            .reorder(c, x, y)
            .vectorize(x, 2*vec, TailStrategy::RoundUp)
            .unroll(c);

        demosaiced->intermed_compute_at.set({processed, yi});
        demosaiced->intermed_store_at.set({processed, yo});
        //demosaiced->output_compute_at.set({curved, x});

        // We can generate slightly better code if we know the splits divide the extent.
        processed
            .bound(c, 0, 3)
            .bound(x, 0, ((out_width)/(2*vec))*(2*vec))
            .bound(y, 0, (out_height/strip_size)*strip_size);
    }
}

HALIDE_REGISTER_GENERATOR(CameraPipe, generator)

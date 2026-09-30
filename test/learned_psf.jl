using SMLMAnalysis
using MicroscopePSFs: MicroscopePSFs
using Test

# See test/long/learned_psf.jl: MicroscopePSFs.save_psf can write a module-qualified type
# tag that load_psf cannot read back, so rewrite it as PSFLearning.save_psf does.
function _save_psf_file(path, psf)
    MicroscopePSFs.save_psf(path, psf)
    SMLMAnalysis.HDF5.h5open(path, "r+") do f
        a = SMLMAnalysis.HDF5.attrs(f)
        a["psf_type"] = last(split(a["psf_type"], '.'))
    end
    return path
end

@testset "DetectFitConfig.psf_file validation" begin
    # These checks run before any SplinePSFModel is built, so they also hold where
    # GaussMLE has no SplinePSFModel (Julia 1.10 CI ignores [sources]). The fit
    # itself is tested in test/long/learned_psf.jl.
    cam = IdealCamera(32, 32, 0.1)
    @test DetectFitConfig().psf_file == ""
    cfg = DetectFitConfig(camera = cam)
    @test SMLMAnalysis._resolve_psf(cfg, cam) === cfg   # no file: unchanged

    mktempdir() do dir
        # A file and a non-default psf_model compete: rejected before the file is read.
        both = DetectFitConfig(
            camera = cam, psf_file = joinpath(dir, "never_read.h5"),
            fitter = GaussMLEConfig(psf_model = GaussianXYNB(0.13f0))
        )
        @test_throws ArgumentError SMLMAnalysis._resolve_psf(both, cam)

        missing_file = DetectFitConfig(camera = cam, psf_file = joinpath(dir, "missing.h5"))
        @test_throws ArgumentError SMLMAnalysis._resolve_psf(missing_file, cam)

        # A file that holds anything but a 3D SplinePSF is rejected, naming the type.
        p2 = _save_psf_file(joinpath(dir, "gauss.h5"), MicroscopePSFs.GaussianPSF(0.13))
        err = try
            SMLMAnalysis._resolve_psf(DetectFitConfig(camera = cam, psf_file = p2), cam)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("GaussianPSF", err.msg)
    end
end

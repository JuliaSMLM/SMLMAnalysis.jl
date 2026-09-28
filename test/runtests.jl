# Lab test runner (admiral decisions 0008, 0009). Identical in every package: do not edit;
# declare groups in test_groups.toml and put tests in the group folders instead.
# GROUP: unset or "Core" runs Core; "QA" or "GPU,Long" runs those; "Everything" runs every
# declared group whose requirement is met on this machine. LAB_TEST_SUMMARY=<path> writes a
# TOML summary there.
using Test, TOML, Pkg

const GROUPS = ("Core", "QA", "GPU", "Data", "Long", "Hardware")
const TESTDIR = @__DIR__
const PKGROOT = dirname(TESTDIR)
const CFG = TOML.parsefile(joinpath(TESTDIR, "test_groups.toml"))

folder(g) = g == "Core" ? TESTDIR : joinpath(TESTDIR, lowercase(g))
testfiles(g) = sort(
    [
        joinpath(folder(g), f) for f in readdir(folder(g))
            if endswith(f, ".jl") && f != "runtests.jl" && isfile(joinpath(folder(g), f))
    ]
)

# Layout check: a test file must never silently not run.
for g in keys(CFG)
    g in GROUPS ||
        error("test_groups.toml: unknown group $g (allowed: $(join(GROUPS, ", ")))")
    isdir(folder(g)) && !isempty(testfiles(g)) ||
        error("group $g is declared but $(folder(g)) has no .jl files")
end
for d in readdir(TESTDIR)
    isdir(joinpath(TESTDIR, d)) || continue
    any(g -> g != "Core" && haskey(CFG, g) && lowercase(g) == d, GROUPS) ||
        error("test/$d/ is not the folder of a group declared in test_groups.toml")
end

# Returns nothing when every requirement of group g is met, else the reason it is not.
function unmet(g)
    for r in get(CFG[g], "requires", String[])
        if r == "cuda"
            Base.find_package("CUDA") === nothing &&
                return "CUDA.jl is not in the $g environment"
            Core.eval(Main, :(import CUDA))
            Base.invokelatest(() -> Main.CUDA.functional()) ||
                return "CUDA.functional() is false"
        elseif r == "data"
            p = expanduser(get(CFG[g], "data_path", ""))
            isempty(p) && error("test_groups.toml: $g requires data but sets no data_path")
            ispath(p) || return "data_path $p does not exist"
        elseif r == "hardware"
            get(ENV, "TEST_HARDWARE", "") == "1" || return "TEST_HARDWARE=1 is not set"
        else
            error(
                "test_groups.toml: $g has unknown requirement \"$r\" (cuda, data, hardware)"
            )
        end
    end
    return nothing
end

# A group with its own test/<group>/Project.toml runs in a temporary copy of that
# environment with the package developed into it, so its dependencies (CUDA) never
# enter the Core env.
function with_group_env(f, g)
    proj = joinpath(folder(g), "Project.toml")
    (g != "Core" && isfile(proj)) || return f()
    prev, tmp = Base.active_project(), mktempdir()
    cp(proj, joinpath(tmp, "Project.toml"))
    Pkg.activate(tmp; io = devnull)
    try
        Pkg.develop(Pkg.PackageSpec(path = PKGROOT); io = devnull)
        Pkg.instantiate(; io = devnull)
        return f()
    finally
        Pkg.activate(prev; io = devnull)
    end
end

function counts(ts)
    c = Test.get_test_counts(ts)  # a Tuple before Julia 1.11, a TestCounts after
    return c isa Tuple ? (pass = c[1] + c[5], fail = 0, error = 0, broken = c[4] + c[8]) :
        (
            pass = c.passes + c.cumulative_passes, fail = 0, error = 0,
            broken = c.broken + c.cumulative_broken,
        )
end

# Runs group g: each file in its own module (like SafeTestsets) and its own @testset.
function rungroup(g)
    t0 = time()
    return with_group_env(g) do
        reason = unmet(g)
        reason === nothing ||
            return Dict{String, Any}("ran" => false, "passed" => false, "reason" => reason)
        printstyled("GROUP $g\n"; bold = true)
        c = try
            counts(
                @testset "$g" begin
                    for f in testfiles(g)
                        @testset "$(relpath(f, TESTDIR))" begin
                            Core.eval(
                                Main, :(
                                    module $(gensym(:testfile))
                                    include($f)
                                    end
                                )
                            )
                        end
                    end
                end
            )
        catch e
            e isa Test.TestSetException || rethrow()
            (pass = e.pass, fail = e.fail, error = e.error, broken = e.broken)
        end
        Dict{String, Any}(
            "ran" => true, "passed" => c.fail + c.error == 0, "pass" => c.pass,
            "fail" => c.fail, "error" => c.error, "broken" => c.broken,
            "seconds" => round(time() - t0; digits = 1)
        )
    end
end

sel = strip(get(ENV, "GROUP", ""))
sel = isempty(sel) ? "Core" : sel
explicit = sel != "Everything"
wanted = explicit ? strip.(split(sel, ",")) : [g for g in GROUPS if haskey(CFG, g)]
for g in wanted
    haskey(CFG, g) || error("GROUP=$sel: $g is not declared in test_groups.toml")
end

results = Dict{String, Any}()
for g in GROUPS
    g in wanted || continue
    results[g] = rungroup(g)
    r = results[g]
    r["ran"] || println(explicit ? "ERROR" : "SKIPPED", " group $g: ", r["reason"])
end

if haskey(ENV, "LAB_TEST_SUMMARY")
    open(ENV["LAB_TEST_SUMMARY"], "w") do io
        TOML.print(
            io, Dict(
                "julia" => string(VERSION), "host" => first(split(gethostname(), '.')),
                "selection" => sel, "groups" => results
            ); sorted = true
        )
    end
end
bad = [g for (g, r) in results if r["ran"] ? !r["passed"] : explicit]
isempty(bad) || error("test groups failed or could not run: $(join(sort(bad), ", "))")

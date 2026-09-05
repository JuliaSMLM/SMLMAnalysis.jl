# ============================================================
# AI coding-assistant guide installer (lab skills-installer convention)
# ============================================================
# `install_agent_guide()` writes a hierarchical, version-stamped guide to the whole
# JuliaSMLM ecosystem into a repo (or the user's home) so a coding assistant (Claude
# Code or Codex) has the APIs of SMLMAnalysis and every sub-package on hand.
#
# Content is assembled at call time from the versions RESOLVED in the current
# environment — each package's `api_overview.md` (README fallback) is read from its
# `pkgdir`, so the guide never drifts from what is actually installed.
#
# Follows the lab-wide convention for package-shipped skill/agent-guide installers:
#   - namespaced install dir  .claude/skills/<pkgprefix>-<skill>/  (collision-safe)
#   - provenance stamp (x-installer / x-source-version / x-source-commit /
#     x-installed-format) written into the installed SKILL.md frontmatter
#   - track=false gitignore default; copy-never-symlink
#   - re-run refreshes only THIS installer's own stamped install; a foreign or
#     hand-made target is refused unless overwrite=true
#   - stamp-scoped uninstall (uninstall_agent_guide) + a doctor (agent_guide_status)
#   - Codex: one managed, stamped block per installing package in AGENTS.md

const _INSTALLER        = "SMLMAnalysis"
const _SKILL_DIRNAME    = "smlma-ecosystem"   # <pkgprefix>-<skill>
const _INSTALLED_FORMAT = 1
const _CODEX_BUNDLE     = "smlm-agent-guide"

# Ordered ecosystem package list (module => one-line role). SMLMAnalysis itself is
# prepended in _collect_pkg_docs as the integration layer. Order mirrors the
# dependency hierarchy: core types first, then single-step packages, then the
# higher-level grouping/clustering packages.
_ecosystem_packages() = (
    (SMLMData,            "Core types shared by every package: Emitter2D/3D(Fit), Camera, BasicSMLD, ROIBatch"),
    (SMLMBoxer,           "ROI detection from raw camera images (difference-of-Gaussians)"),
    (GaussMLE,            "GPU-accelerated maximum-likelihood PSF fitting (Gaussian/astigmatic models)"),
    (MicroscopePSFs,      "PSF models: Gaussian, Airy, spline, and vector PSFs"),
    (SMLMSim,             "SMLM data simulation and fluorophore blinking kinetics"),
    (SMLMFrameConnection, "Linking localizations across frames + uncertainty calibration"),
    (SMLMDriftCorrection, "Sample-drift correction (entropy-based) and cross-channel alignment"),
    (SMLMRender,          "Super-resolution image rendering (histogram/Gaussian/circle/ellipse)"),
    (SMLMClustering,      "Clustering backends, spatial-tendency statistics, and edge classification"),
    (SMLMBaGoL,           "Bayesian grouping of localizations into high-precision emitters (RJMCMC)"),
)

# Faithful to the hierarchy in CLAUDE.md / README.
_dependency_tree() = """
SMLMData  (core types — no deps)
   |
   +-- SMLMBoxer ----+
   +-- GaussMLE -----+--- detection + fitting
   +-- MicroscopePSFs --- SMLMSim   (simulation)
   +-- SMLMFrameConnection
   +-- SMLMDriftCorrection
   +-- SMLMRender
   +-- SMLMClustering ---- SMLMBaGoL (Bayesian grouping)
   |
SMLMAnalysis  (integrates all of the above into one analyze() pipeline)
"""

# Single-line YAML description for the Claude skill frontmatter (also the lead line of
# the Codex block). Kept free of `: ` so it is a safe YAML plain scalar.
_guide_description() = "JuliaSMLM SMLM analysis ecosystem — the unified analyze() pipeline (detect, fit, filter, frame-connect, drift-correct, render, cluster, BaGoL) plus the full APIs of every sub-package. Use when writing or debugging single-molecule-localization-microscopy code with SMLMAnalysis.jl or any JuliaSMLM package."

# --- provenance stamp ----------------------------------------------------------

_source_version() = try string(pkgversion(@__MODULE__)) catch; "unknown" end

# Short commit of the installing package's checkout, or "none" for a registered
# install with no git tree. Never throws — provenance is best-effort.
function _source_commit()
    root = pkgdir(@__MODULE__)
    root === nothing && return "none"
    try
        c = strip(readchomp(`git -C $root rev-parse --short=7 HEAD`))
        isempty(c) ? "none" : String(c)
    catch
        "none"
    end
end

_stamp_pairs() = ("x-installer"        => _INSTALLER,
                  "x-source-version"   => _source_version(),
                  "x-source-commit"    => _source_commit(),
                  "x-installed-format" => string(_INSTALLED_FORMAT))

_stamp_inline() = join(("$k: $v" for (k, v) in _stamp_pairs()), " | ")

# Read a leading-frontmatter field (between the first two `---` lines). nothing if
# the file/field is absent — used to identify our own installs by x-installer.
function _frontmatter_field(path::AbstractString, key::AbstractString)
    isfile(path) || return nothing
    seen = 0
    for ln in readlines(path)
        s = strip(ln)
        if s == "---"
            seen += 1
            seen == 2 && break
            continue
        end
        if seen == 1 && startswith(s, key * ":")
            return strip(s[length(key)+2:end], ['"', ' '])
        end
    end
    nothing
end

# Read an inline `key: value` stamp field from the Codex GUIDE.md comment header.
function _guide_field(guide::AbstractString, key::AbstractString)
    isfile(guide) || return nothing
    m = match(Regex(key * raw":\s*([^\s|]+)"), read(guide, String))
    m === nothing ? nothing : String(m.captures[1])
end

# --- content collection --------------------------------------------------------

function _pkg_entry(mod::Module, role::AbstractString)
    name = string(nameof(mod))
    root = pkgdir(mod)
    ver  = try string(pkgversion(mod)) catch; "unknown" end
    (text, source) = _read_api_text(root, name)
    (; name, role, version = ver, text, source)
end

# api_overview.md is the ecosystem's AI-parseable reference convention; README is the
# fallback. Returns (text, source_label).
function _read_api_text(root::Union{AbstractString,Nothing}, name::AbstractString)
    root === nothing && return ("_(package directory for $name not found)_\n", "none")
    for fname in ("api_overview.md", "README.md")
        p = joinpath(root, fname)
        isfile(p) && return (read(p, String), fname)
    end
    ("_(no api_overview.md or README.md found for $name)_\n", "none")
end

function _collect_pkg_docs()
    entries = Any[_pkg_entry(@__MODULE__,
        "Integration layer — the unified analyze(state, config) pipeline that orchestrates every package below")]
    for (mod, role) in _ecosystem_packages()
        push!(entries, _pkg_entry(mod, role))
    end
    entries
end

# --- rendering -----------------------------------------------------------------

function _render_guide(entries, refprefix::AbstractString)
    io = IOBuffer()
    println(io, "# SMLMAnalysis — JuliaSMLM ecosystem API guide")
    println(io)
    println(io, "Generated by `SMLMAnalysis.install_agent_guide()` from the package versions ",
                "resolved in this environment. Re-run it to refresh.")
    println(io)
    println(io, "SMLMAnalysis integrates the JuliaSMLM single-molecule-localization-microscopy ",
                "packages into one `analyze(state, config)` pipeline: the config type selects the ",
                "operation via multiple dispatch, and steps compose in any order after detection/fitting. ",
                "Coordinates are in microns throughout; raw data is `(height, width, frames)` image stacks.")
    println(io)
    println(io, "## Dependency hierarchy")
    println(io, "```")
    print(io, _dependency_tree())
    println(io, "```")
    println(io)
    println(io, "## Packages (versions resolved in this environment)")
    println(io)
    for e in entries
        println(io, "- **$(e.name)** `v$(e.version)` — $(e.role)  ")
        println(io, "  Full API: [`$(refprefix)/$(e.name).md`]($(refprefix)/$(e.name).md) *(source: $(e.source))*")
    end
    println(io)
    println(io, "## Where to start")
    println(io, "- Entry point: `analyze(image_stacks, AnalysisConfig(...))`, or step-by-step ",
                "`analyze(state, SomeConfig())`.")
    println(io, "- Full SMLMAnalysis API (pipeline, step configs, I/O, multi-target): see ",
                "[`$(refprefix)/SMLMAnalysis.md`]($(refprefix)/SMLMAnalysis.md).")
    println(io, "- Each `reference/*.md` above is that package's own `api_overview.md`, copied verbatim.")
    String(take!(io))
end

# SKILL.md = stamped YAML frontmatter + the rendered guide body.
function _skill_text(entries)
    io = IOBuffer()
    println(io, "---")
    println(io, "name: ", _SKILL_DIRNAME)
    println(io, "description: \"", _guide_description(), "\"")
    for (k, v) in _stamp_pairs()
        println(io, k, ": ", v)
    end
    println(io, "---")
    println(io)
    print(io, _render_guide(entries, "reference"))
    String(take!(io))
end

function _write_reference!(refdir::AbstractString, e)
    mkpath(refdir)
    p = joinpath(refdir, e.name * ".md")
    open(p, "w") do io
        println(io, "<!-- $(e.name) v$(e.version) — copied from $(e.source) by ",
                    "SMLMAnalysis.install_agent_guide(). Do not edit; re-run to refresh. -->")
        println(io)
        print(io, e.text)
    end
    p
end

# Remove only the reference files this installer generated (per the manifest written at
# install time, falling back to the known package names for a pre-manifest install),
# then drop reference/ only if nothing else is in it. Never follows symlinks.
function _remove_generated_reference!(refdir::AbstractString)
    isdir(refdir) || return nothing
    islink(refdir) && return nothing
    manifest = joinpath(refdir, ".smlma-manifest")
    names = if isfile(manifest)
        filter(!isempty, strip.(readlines(manifest)))
    else
        vcat([string(_INSTALLER) * ".md"], [string(nameof(mod)) * ".md" for (mod, _) in _ecosystem_packages()])
    end
    for n in names
        p = joinpath(refdir, n)
        isfile(p) && !islink(p) && rm(p)
    end
    isfile(manifest) && rm(manifest)
    if isempty(readdir(refdir))
        rm(refdir)
    else
        @info "$(refdir): kept files this installer did not write" remaining=readdir(refdir)
    end
    nothing
end

# --- git / AGENTS.md side effects ----------------------------------------------

# Append patterns to <root>/.gitignore (creating it if needed), skipping any already
# present. Preserves a missing trailing newline in the existing file.
function _ensure_gitignored!(root::AbstractString, patterns)
    gi   = joinpath(root, ".gitignore")
    raw  = isfile(gi) ? read(gi, String) : ""
    have = Set(strip.(split(raw, '\n')))
    todo = [p for p in patterns if !(strip(p) in have)]
    isempty(todo) && return gi
    buf = raw
    isempty(buf) || endswith(buf, '\n') || (buf *= '\n')
    isempty(buf) || (buf *= '\n')
    buf *= "# SMLMAnalysis agent guide — installed with track=false (install_agent_guide)\n"
    buf *= join(todo, '\n') * '\n'
    write(gi, buf)
    gi
end

# Package-scoped markers so multiple packages can each own a block in one AGENTS.md.
const _AGENTS_BEGIN = "<!-- BEGIN SMLMAnalysis agent-guide (managed by install_agent_guide) -->"
const _AGENTS_END   = "<!-- END SMLMAnalysis agent-guide -->"

# Insert or replace our delimited block in AGENTS.md, leaving other content untouched.
# Idempotent: re-running refreshes the block in place.
function _upsert_agents_block!(path::AbstractString, block::AbstractString)
    raw     = isfile(path) ? read(path, String) : ""
    managed = string(_AGENTS_BEGIN, '\n', block, '\n', _AGENTS_END)
    b = findfirst(_AGENTS_BEGIN, raw)
    e = findfirst(_AGENTS_END, raw)
    if b !== nothing && e !== nothing && first(e) > first(b)
        newraw = raw[1:prevind(raw, first(b))] * managed * raw[nextind(raw, last(e)):end]
        write(path, newraw)
    else
        buf = raw
        isempty(buf) || endswith(buf, '\n') || (buf *= '\n')
        isempty(buf) || (buf *= '\n')
        write(path, buf * managed * '\n')
    end
    path
end

# Remove our block from AGENTS.md; returns true if a block was present and removed.
#
# Exact inverse of _upsert_agents_block!, not a generic strip: that function, on first
# insert, normalizes any nonempty pre-existing content to end with a blank line (one
# newline to terminate its last line if needed, plus one blank-line newline) before the
# BEGIN marker, and always places exactly one literal "\n" right after the END marker
# (a later refresh copies both sides verbatim, so this holds no matter how many times
# the block has been refreshed). Removal undoes only those two installer-added
# characters — one trailing "\n" off `pre`, one leading "\n" off `post` — and returns
# `pre * post` otherwise byte-for-byte, so user content (including indentation the
# user appended after the block) is never touched.
function _remove_agents_block!(path::AbstractString)
    isfile(path) || return false
    raw = read(path, String)
    b = findfirst(_AGENTS_BEGIN, raw)
    e = findfirst(_AGENTS_END, raw)
    (b === nothing || e === nothing || first(e) < first(b)) && return false
    pre  = raw[1:prevind(raw, first(b))]
    post = raw[nextind(raw, last(e)):end]
    endswith(pre, "\n\n") && (pre = chop(pre))
    startswith(post, "\n") && (post = chop(post; head = 1, tail = 0))
    write(path, pre * post)
    true
end

_codex_block() = string(
    "<!-- ", _stamp_inline(), " -->\n",
    _guide_description(), "\n\n",
    "When working with JuliaSMLM / SMLMAnalysis code, read the ecosystem guide at ",
    "`", _CODEX_BUNDLE, "/GUIDE.md` and the per-package API references under ",
    "`", _CODEX_BUNDLE, "/reference/`.")

# expanduser is a no-op on Windows; handle a leading ~ ourselves on every platform.
function _expand_home(dir::AbstractString)
    d = String(dir)
    d == "~" && return homedir()
    (startswith(d, "~/") || startswith(d, "~\\")) && return joinpath(homedir(), d[3:end])
    return expanduser(d)
end

# Skill/bundle directory for a (tool, scope, dir) triple.
function _install_dir(tool::Symbol, scope::Symbol, dir::AbstractString)
    if tool == :claude
        scope == :project ? joinpath(dir, ".claude", "skills", _SKILL_DIRNAME) :
                            joinpath(homedir(), ".claude", "skills", _SKILL_DIRNAME)
    else
        base = scope == :project ? String(dir) : joinpath(homedir(), ".codex")
        joinpath(base, _CODEX_BUNDLE)
    end
end

# --- public entry points -------------------------------------------------------

"""
    install_agent_guide(; tool=:claude, scope=:project, track=false,
                          overwrite=false, dir=pwd()) -> String

Install a hierarchical, version-stamped guide to the JuliaSMLM ecosystem for an AI
coding assistant, so it has the APIs of `SMLMAnalysis` and every sub-package on hand.

The guide is assembled at call time from the versions resolved in the current
environment — each package's `api_overview.md` (falling back to `README.md`) is read
from its `pkgdir` and copied into a `reference/` bundle, with a top-level map linking
them. Re-running refreshes it against whatever is currently installed.

Follows the lab skills-installer convention: the install carries a provenance stamp
(`x-installer`, `x-source-version`, `x-source-commit`), so re-running refreshes only
*this* installer's own install, and [`uninstall_agent_guide`](@ref) /
[`agent_guide_status`](@ref) act only on stamped installs. A refresh regenerates only
the reference files listed in its own manifest (`reference/.smlma-manifest`), so files
you add inside `reference/` yourself survive; a symlinked target, wrapper, or reference
directory is refused with an `ArgumentError` rather than mutated through.

# Keyword arguments
- `tool::Symbol = :claude` — target assistant:
  - `:claude` writes a Claude Code **skill** (`.claude/skills/smlma-ecosystem/SKILL.md`
    + `reference/*.md`).
  - `:codex` writes a `smlm-agent-guide/` bundle and a managed block in `AGENTS.md`.
- `scope::Symbol = :project` — `:project` → into `dir` (the repo); `:user` → into your
  home (`~/.claude` / `~/.codex`, applying to every project you open).
- `track::Bool = false` — **project scope only.** When `false` (default) the installed
  files are added to the repo's `.gitignore`, keeping the guide out of history. Pass
  `track=true` to commit and share it. Ignored (with a warning) for `scope=:user`.
- `overwrite::Bool = false` — replace a target that was **not** installed by this
  installer (a hand-made skill, or another package's install). Refreshing our own
  stamped install never needs it.
- `dir::AbstractString = pwd()` — the repo root for `scope=:project`.

Returns the path of the installed skill/bundle directory.

# Examples
```julia
using SMLMAnalysis
install_agent_guide()                    # Claude skill in ./.claude, gitignored
install_agent_guide(track=true)          # …and committed to the repo
install_agent_guide(tool=:codex)         # Codex AGENTS.md + bundle in this repo
install_agent_guide(scope=:user)         # Claude skill for all your projects
```
"""
function install_agent_guide(; tool::Symbol = :claude,
                               scope::Symbol = :project,
                               track::Bool = false,
                               overwrite::Bool = false,
                               dir::AbstractString = pwd())
    tool in (:claude, :codex) ||
        throw(ArgumentError("tool must be :claude or :codex, got :$tool"))
    scope in (:project, :user) ||
        throw(ArgumentError("scope must be :project or :user, got :$scope"))
    if track && scope == :user
        @warn "track applies only to scope=:project (the :user guide lives outside any repo); ignoring track=true"
    end
    dir = _expand_home(dir)   # `dir="~/repo"` must not create a literal `~` under the cwd

    entries   = _collect_pkg_docs()
    gitignore = !track && scope == :project
    target    = _install_dir(tool, scope, dir)

    if tool == :claude
        wrapper = joinpath(target, "SKILL.md")
        refdir  = joinpath(target, "reference")
        # Never mutate through a symlink: a symlinked target/wrapper/reference dir may
        # point outside the install dir entirely (uninstall's old rm(target) only ever
        # unlinked the link, but write(wrapper, ...) and rm(refdir; recursive=true)
        # follow a symlink and can clobber or delete files that live elsewhere).
        islink(target) && throw(ArgumentError("$target is a symlink; refusing to install through it"))
        islink(wrapper) && throw(ArgumentError("$wrapper is a symlink; refusing to install through it"))
        islink(refdir) && throw(ArgumentError("$refdir is a symlink; refusing to install through it"))
        # Own-install idempotency: refresh ours freely; refuse a foreign/unstamped
        # target unless overwrite=true. Guard whenever the dir exists with content —
        # NOT only when SKILL.md is present — so a hand-made dir with a reference/ but
        # no SKILL.md is not silently wiped by the rm below.
        if isdir(target) && !isempty(readdir(target))
            owner = isfile(wrapper) ? _frontmatter_field(wrapper, "x-installer") : nothing
            owner == _INSTALLER || overwrite ||
                error("$target already exists and was not installed by $_INSTALLER " *
                      "(x-installer=$(owner === nothing ? "none" : owner)). Pass overwrite=true to replace it.")
        end
        # Regenerate only the reference files we generated last time (per the
        # manifest); anything a user added inside reference/ survives a refresh.
        _remove_generated_reference!(refdir)
        mkpath(refdir)
        for e in entries
            _write_reference!(refdir, e)
        end
        write(joinpath(refdir, ".smlma-manifest"), join((e.name * ".md" for e in entries), '\n') * '\n')
        write(wrapper, _skill_text(entries))
        gitignore && _ensure_gitignored!(dir, ["/.claude/skills/$_SKILL_DIRNAME/"])
        return target
    else # :codex
        guide  = joinpath(target, "GUIDE.md")
        refdir = joinpath(target, "reference")
        islink(target) && throw(ArgumentError("$target is a symlink; refusing to install through it"))
        islink(guide) && throw(ArgumentError("$guide is a symlink; refusing to install through it"))
        islink(refdir) && throw(ArgumentError("$refdir is a symlink; refusing to install through it"))
        # Same guard as :claude — refuse any non-empty foreign target, not only one
        # that already has our GUIDE.md, so a stray reference/ is not wiped.
        if isdir(target) && !isempty(readdir(target))
            owner = isfile(guide) ? _guide_field(guide, "x-installer") : nothing
            owner == _INSTALLER || overwrite ||
                error("$target already exists and was not installed by $_INSTALLER. Pass overwrite=true to replace it.")
        end
        _remove_generated_reference!(refdir)
        mkpath(refdir)
        for e in entries
            _write_reference!(refdir, e)
        end
        write(joinpath(refdir, ".smlma-manifest"), join((e.name * ".md" for e in entries), '\n') * '\n')
        write(guide, string("<!-- ", _stamp_inline(), " -->\n\n", _render_guide(entries, "reference")))
        agents = joinpath(dirname(target), "AGENTS.md")
        _upsert_agents_block!(agents, _codex_block())
        gitignore && _ensure_gitignored!(dir, ["/$_CODEX_BUNDLE/"])
        return target
    end
end

"""
    uninstall_agent_guide(; tool=:claude, scope=:project, dir=pwd()) -> Vector{String}

Remove a guide previously installed by SMLMAnalysis. Acts **only** on a target carrying
SMLMAnalysis's own provenance stamp — a hand-made skill or another package's install is
left untouched — and even then removes only `SKILL.md`/`GUIDE.md`, the reference files
listed in `reference/.smlma-manifest` (falling back to the known package names for a
pre-manifest install), and — for `tool=:codex` — this package's managed block in
`AGENTS.md`. Anything else, including files you added inside `reference/` yourself, is
preserved; `reference/` and the install directory are removed only when nothing else
remains in them, and a non-empty directory is left in place and reported with `@info`.
A symlinked target (or a symlinked `SKILL.md`/`GUIDE.md`) is never followed — it is left
untouched with a `@warn`. Returns the paths removed (empty if nothing of ours was
found, including when the target is a symlink).
"""
function uninstall_agent_guide(; tool::Symbol = :claude,
                                 scope::Symbol = :project,
                                 dir::AbstractString = pwd())
    tool in (:claude, :codex) ||
        throw(ArgumentError("tool must be :claude or :codex, got :$tool"))
    scope in (:project, :user) ||
        throw(ArgumentError("scope must be :project or :user, got :$scope"))
    dir = _expand_home(dir)

    removed = String[]
    target  = _install_dir(tool, scope, dir)
    if tool == :claude
        stamp = joinpath(target, "SKILL.md")
        if islink(target) || islink(stamp)
            @warn "$target is a symlink (or its SKILL.md is); leaving it untouched"
        elseif isdir(target) && _frontmatter_field(stamp, "x-installer") == _INSTALLER
            push!(removed, _remove_own_files!(target, stamp))
        end
    else
        stamp = joinpath(target, "GUIDE.md")
        if islink(target) || islink(stamp)
            @warn "$target is a symlink (or its GUIDE.md is); leaving it untouched"
        elseif isdir(target) && _guide_field(stamp, "x-installer") == _INSTALLER
            push!(removed, _remove_own_files!(target, stamp))
        end
        agents = joinpath(dirname(target), "AGENTS.md")
        _remove_agents_block!(agents) && push!(removed, agents)
    end
    removed
end

# Remove exactly what install_agent_guide writes into `target` — the stamped wrapper
# file and the reference files listed in its manifest — then drop the directory only
# if it is empty. Never `rm(target; recursive=true)`: after an `overwrite=true` install
# into a hand-made directory the stamp is ours but the directory may still hold the
# user's own files, and a recursive delete would take them with it.
function _remove_own_files!(target::AbstractString, stamp::AbstractString)
    isfile(stamp) && rm(stamp)
    refdir = joinpath(target, "reference")
    _remove_generated_reference!(refdir)
    if isempty(readdir(target))
        rm(target)
        return String(target)
    end
    @info "uninstall_agent_guide: removed SMLMAnalysis's files from $target but kept the directory — it still holds files this installer did not write." remaining=readdir(target)
    return String(target)
end

"""
    agent_guide_status(; tool=:claude, scope=:project, dir=pwd()) -> NamedTuple

Doctor: report an installed guide's state — `(; installed, path, source_version,
source_commit, current_version, stale)`. `stale` is true when the stamped
`source_version` differs from the version currently resolved in this environment
(re-run [`install_agent_guide`](@ref) to refresh).
"""
function agent_guide_status(; tool::Symbol = :claude,
                              scope::Symbol = :project,
                              dir::AbstractString = pwd())
    tool in (:claude, :codex) ||
        throw(ArgumentError("tool must be :claude or :codex, got :$tool"))
    scope in (:project, :user) ||
        throw(ArgumentError("scope must be :project or :user, got :$scope"))

    current = _source_version()
    target  = _install_dir(tool, scope, _expand_home(dir))
    stampfile = tool == :claude ? joinpath(target, "SKILL.md") : joinpath(target, "GUIDE.md")
    reader    = tool == :claude ? _frontmatter_field : _guide_field

    # A symlinked target/stamp is never followed elsewhere in this file; report it as
    # not-installed here too rather than reading through it.
    installed = !islink(target) && !islink(stampfile) && reader(stampfile, "x-installer") == _INSTALLER
    sv = installed ? reader(stampfile, "x-source-version") : nothing
    sc = installed ? reader(stampfile, "x-source-commit") : nothing
    (; installed,
       path = target,
       source_version = sv,
       source_commit = sc,
       current_version = current,
       stale = installed && sv !== nothing && sv != current)
end

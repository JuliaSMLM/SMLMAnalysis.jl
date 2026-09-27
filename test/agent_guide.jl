using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "install_agent_guide" begin
    # Arg validation
    @test_throws ArgumentError install_agent_guide(tool = :bogus)
    @test_throws ArgumentError install_agent_guide(scope = :bogus)
    @test_throws ArgumentError uninstall_agent_guide(tool = :bogus)
    @test_throws ArgumentError agent_guide_status(scope = :bogus)

    # Claude, project scope, default track=false → gitignored, namespaced, stamped.
    mktempdir() do dir
        skill = install_agent_guide(dir = dir)
        @test skill == joinpath(dir, ".claude", "skills", "smlma-ecosystem")
        @test isfile(joinpath(skill, "SKILL.md"))

        refs = readdir(joinpath(skill, "reference"))
        # SMLMAnalysis + the 10 ecosystem packages, each with a reference file.
        @test length(refs) == 11
        for name in ("SMLMAnalysis", "SMLMData", "GaussMLE", "SMLMBaGoL", "SMLMRender")
            @test "$name.md" in refs
        end

        skilltext = read(joinpath(skill, "SKILL.md"), String)
        @test occursin("name: smlma-ecosystem", skilltext)
        @test occursin("description:", skilltext)
        @test occursin("x-installer: SMLMAnalysis", skilltext)   # provenance stamp
        @test occursin("x-source-version:", skilltext)
        @test occursin("Dependency hierarchy", skilltext)

        # track=false (default) gitignores the namespaced skill dir.
        @test occursin(".claude/skills/smlma-ecosystem/", read(joinpath(dir, ".gitignore"), String))

        # The rewritten .gitignore follows the umask, like a direct `write`
        # would — not the fixed 0600 a mktemp-created temp carries over
        # (the bug _replace_atomically's private-temp-dir + write() path avoids).
        control = joinpath(dir, "control.txt")
        write(control, "control")
        @test filemode(joinpath(dir, ".gitignore")) & 0o777 == filemode(control) & 0o777

        # Doctor: freshly installed, not stale.
        st = agent_guide_status(dir = dir)
        @test st.installed
        @test !st.stale
        @test st.source_version == st.current_version

        # Own-install idempotency: re-running our OWN stamped install refreshes
        # WITHOUT overwrite (no error), returning the same path.
        @test install_agent_guide(dir = dir) == skill

        # A foreign/unstamped skill in the same dir IS refused unless overwrite.
        write(joinpath(skill, "SKILL.md"), "---\nname: smlma-ecosystem\n---\nhand-made\n")
        @test_throws ErrorException install_agent_guide(dir = dir)
        @test install_agent_guide(dir = dir, overwrite = true) == skill
        @test occursin("x-installer: SMLMAnalysis", read(joinpath(skill, "SKILL.md"), String))

        # Uninstall removes only our stamped install.
        @test uninstall_agent_guide(dir = dir) == [skill]
        @test !isdir(skill)
        @test !agent_guide_status(dir = dir).installed
        @test isempty(uninstall_agent_guide(dir = dir))   # nothing left → no-op
    end

    # Guard-gap regression: a hand-made target dir with a reference/ but NO
    # SKILL.md must be REFUSED (not silently wiped) unless overwrite=true.
    mktempdir() do dir
        skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
        mkpath(joinpath(skill, "reference"))
        write(joinpath(skill, "reference", "keep.md"), "hand-made")
        @test_throws ErrorException install_agent_guide(dir = dir)
        @test isfile(joinpath(skill, "reference", "keep.md"))            # survived
        @test install_agent_guide(dir = dir, overwrite = true) == skill  # overwrite proceeds
        @test isfile(joinpath(skill, "SKILL.md"))
    end

    # Uninstall leaves a foreign skill untouched.
    mktempdir() do dir
        skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
        mkpath(skill)
        write(joinpath(skill, "SKILL.md"), "---\nname: smlma-ecosystem\n---\nnot ours\n")
        @test isempty(uninstall_agent_guide(dir = dir))
        @test isdir(skill)
    end

    # Regression (data loss): overwrite=true onto a hand-made dir stamps it as ours
    # while preserving the user's own files — a later uninstall must remove ONLY
    # what we wrote (SKILL.md + reference/) and leave the still-populated directory
    # in place, never `rm -r` the user's notes.md / scripts/ along with it.
    mktempdir() do dir
        skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
        mkpath(joinpath(skill, "scripts"))
        write(joinpath(skill, "SKILL.md"), "---\nname: smlma-ecosystem\n---\nhand-made\n")
        write(joinpath(skill, "notes.md"), "my notes")
        write(joinpath(skill, "scripts", "tool.py"), "print(1)")
        @test install_agent_guide(dir = dir, overwrite = true) == skill
        @test isfile(joinpath(skill, "notes.md"))                       # install preserved it
        @test uninstall_agent_guide(dir = dir) == [skill]
        @test isdir(skill)                                              # dir kept (non-empty)
        @test isfile(joinpath(skill, "notes.md"))                       # user's files survive
        @test isfile(joinpath(skill, "scripts", "tool.py"))
        @test !isfile(joinpath(skill, "SKILL.md"))                      # ours removed
        @test !isdir(joinpath(skill, "reference"))
        @test !agent_guide_status(dir = dir).installed
        @test isempty(uninstall_agent_guide(dir = dir))                 # nothing of ours left
    end

    # Same contract for the Codex bundle: user files inside smlm-agent-guide/ survive.
    mktempdir() do dir
        bundle = joinpath(dir, "smlm-agent-guide")
        mkpath(bundle)
        write(joinpath(bundle, "GUIDE.md"), "hand-made\n")
        write(joinpath(bundle, "mine.txt"), "keep")
        @test install_agent_guide(dir = dir, tool = :codex, overwrite = true) == bundle
        removed = uninstall_agent_guide(dir = dir, tool = :codex)
        @test bundle in removed
        @test isdir(bundle) && isfile(joinpath(bundle, "mine.txt"))
        @test !isfile(joinpath(bundle, "GUIDE.md")) && !isdir(joinpath(bundle, "reference"))
    end

    # reference/ survival across a refresh AND an uninstall (Claude): a file the
    # user adds inside reference/ is never touched, but files we generated
    # (identified by their per-file provenance stamp) ARE regenerated by a
    # refresh / removed by uninstall.
    mktempdir() do dir
        skill = install_agent_guide(dir = dir)
        refdir = joinpath(skill, "reference")
        notes = joinpath(refdir, "personal-notes.md")
        write(notes, "mine")

        datamd = joinpath(refdir, "SMLMData.md")
        rm(datamd)   # force a visible change so we can tell a refresh regenerates it

        @test install_agent_guide(dir = dir) == skill      # plain refresh, no overwrite
        @test read(notes, String) == "mine"                # user file survived untouched
        @test isfile(datamd)                               # ours was regenerated

        removed = uninstall_agent_guide(dir = dir)
        @test removed == [skill]
        @test isfile(notes)                                 # user file still there
        @test isdir(refdir)                                 # dir kept (notes still inside)
        @test !isfile(datamd)                               # generated files removed
        @test !isfile(joinpath(skill, "SKILL.md"))
        @test !agent_guide_status(dir = dir).installed
    end

    # Same contract for Codex.
    mktempdir() do dir
        bundle = install_agent_guide(dir = dir, tool = :codex)
        refdir = joinpath(bundle, "reference")
        notes = joinpath(refdir, "personal-notes.md")
        write(notes, "mine")

        datamd = joinpath(refdir, "SMLMData.md")
        rm(datamd)

        @test install_agent_guide(dir = dir, tool = :codex) == bundle
        @test read(notes, String) == "mine"
        @test isfile(datamd)

        removed = uninstall_agent_guide(dir = dir, tool = :codex)
        @test bundle in removed
        @test isfile(notes)
        @test isdir(refdir)
        @test !isfile(datamd)
        @test !isfile(joinpath(bundle, "GUIDE.md"))
        @test !agent_guide_status(dir = dir, tool = :codex).installed
    end

    # Data-loss regression: a DIRECTORY (not just a symlink) sitting at a reference
    # file's name must be refused, never wiped by mv(...; force=true)'s implicit
    # recursive rm of the destination. Exercised via :claude; :codex shares the same
    # _preflight_reference_files helper.
    mktempdir() do dir
        skill = install_agent_guide(dir = dir)
        refdir = joinpath(skill, "reference")
        datamd = joinpath(refdir, "SMLMData.md")
        rm(datamd)
        mkpath(datamd)
        write(joinpath(datamd, "keep.txt"), "keep")

        @test_throws ArgumentError install_agent_guide(dir = dir)
        @test isdir(datamd)                                 # not wiped
        @test isfile(joinpath(datamd, "keep.txt"))           # contents survived
    end

    # Same regression for AGENTS.md itself: a directory there must be refused
    # BEFORE any bundle content is written (nothing half-installed).
    mktempdir() do dir
        mkpath(joinpath(dir, "AGENTS.md"))
        @test_throws ArgumentError install_agent_guide(dir = dir, tool = :codex)
        @test !ispath(joinpath(dir, "smlm-agent-guide"))    # bundle never written
    end

    # agent_guide's writes now go through _replace_atomically (src/io/atomic.jl),
    # same as save_smld: a failed write (destination is a directory, refused
    # outright) must leave no stray temp file behind.
    mktempdir() do dir
        mkpath(joinpath(dir, "somedir"))
        @test_throws ArgumentError SMLMAnalysis._replace_atomically(joinpath(dir, "somedir")) do tmp
            write(tmp, "x")
        end
        @test readdir(dir) == ["somedir"]   # no stray temp file
    end

    if Sys.isunix()
        # (a) target itself is a symlink to a real, stamped install elsewhere:
        # neither install nor uninstall may follow it (Claude).
        mktempdir() do root
            actual_repo = joinpath(root, "actual_repo")
            mkpath(actual_repo)
            actual_target = install_agent_guide(dir = actual_repo)

            repo = joinpath(root, "repo")
            mkpath(joinpath(repo, ".claude", "skills"))
            target = joinpath(repo, ".claude", "skills", "smlma-ecosystem")
            symlink(actual_target, target; dir_target = true)

            @test isempty(uninstall_agent_guide(dir = repo))
            @test isfile(joinpath(actual_target, "SKILL.md"))   # real install untouched
            @test_throws ArgumentError install_agent_guide(dir = repo)
        end

        # (a) same for Codex.
        mktempdir() do root
            actual_repo = joinpath(root, "actual_repo")
            mkpath(actual_repo)
            actual_bundle = install_agent_guide(dir = actual_repo, tool = :codex)

            repo = joinpath(root, "repo")
            mkpath(repo)
            target = joinpath(repo, "smlm-agent-guide")
            symlink(actual_bundle, target; dir_target = true)

            @test isempty(uninstall_agent_guide(dir = repo, tool = :codex))
            @test isfile(joinpath(actual_bundle, "GUIDE.md"))   # real install untouched
            @test_throws ArgumentError install_agent_guide(dir = repo, tool = :codex)
        end

        # (b) the wrapper itself is a symlink pointing outside the install dir
        # (Claude): install must refuse rather than write through it.
        mktempdir() do dir
            skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
            mkpath(skill)
            external = joinpath(dir, "external.md")
            write(external, "not ours")
            symlink(external, joinpath(skill, "SKILL.md"))

            @test_throws ArgumentError install_agent_guide(dir = dir, overwrite = true)
            @test read(external, String) == "not ours"   # external file untouched
        end

        # (c) a reference file itself is a symlink to an external file: a plain
        # refresh must refuse rather than write through it.
        mktempdir() do dir
            skill = install_agent_guide(dir = dir)
            refdir = joinpath(skill, "reference")
            external = joinpath(dir, "external-ref.md")
            write(external, "not ours")
            rm(joinpath(refdir, "SMLMData.md"))
            symlink(external, joinpath(refdir, "SMLMData.md"))

            @test_throws ArgumentError install_agent_guide(dir = dir)
            @test read(external, String) == "not ours"   # external file untouched
        end

        # (d) the user replaces a generated reference file with their own content
        # (no provenance header): refresh without overwrite refuses and leaves it
        # untouched; overwrite=true replaces it and restores the stamp.
        mktempdir() do dir
            skill = install_agent_guide(dir = dir)
            refdir = joinpath(skill, "reference")
            datamd = joinpath(refdir, "SMLMData.md")
            write(datamd, "hand-edited, no stamp\n")

            @test_throws ErrorException install_agent_guide(dir = dir)
            @test read(datamd, String) == "hand-edited, no stamp\n"

            @test install_agent_guide(dir = dir, overwrite = true) == skill
            @test occursin("SMLMAnalysis.install_agent_guide()", read(datamd, String))
        end

        # (e) a hardlinked wrapper: install must replace it via rename, never
        # truncate it in place, so the OTHER hardlink (the user's) is untouched.
        mktempdir() do dir
            skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
            mkpath(skill)
            external = joinpath(dir, "external-skill.md")
            write(external, "not ours")
            @test ccall(:link, Cint, (Cstring, Cstring), external, joinpath(skill, "SKILL.md")) == 0

            @test install_agent_guide(dir = dir, overwrite = true) == skill
            @test occursin("x-installer: SMLMAnalysis", read(joinpath(skill, "SKILL.md"), String))
            @test read(external, String) == "not ours"   # other hardlink untouched
        end

        # (f) symlinked .gitignore / AGENTS.md: never written through.
        mktempdir() do dir
            external = joinpath(dir, "external-gitignore")
            write(external, "external\n")
            symlink(external, joinpath(dir, ".gitignore"))
            @test_throws ArgumentError install_agent_guide(dir = dir)
            @test read(external, String) == "external\n"
        end
        mktempdir() do dir
            external = joinpath(dir, "external-agents")
            write(external, "external\n")
            symlink(external, joinpath(dir, "AGENTS.md"))
            @test_throws ArgumentError install_agent_guide(dir = dir, tool = :codex)
            @test read(external, String) == "external\n"
        end

        # (g) Codex: the bundle directory is replaced by a symlink to a copy
        # elsewhere — uninstall must not follow it, and must leave the AGENTS.md
        # block in place (it must not go on to strip the block once the bundle
        # itself is judged unsafe).
        mktempdir() do dir
            bundle = install_agent_guide(dir = dir, tool = :codex)
            copy_elsewhere = joinpath(dir, "bundle-copy")
            cp(bundle, copy_elsewhere)
            rm(bundle; recursive = true)
            symlink(copy_elsewhere, bundle; dir_target = true)

            @test isempty(uninstall_agent_guide(dir = dir, tool = :codex))
            agents = read(joinpath(dir, "AGENTS.md"), String)
            @test occursin("BEGIN SMLMAnalysis agent-guide", agents)
        end

        # (h) a symlinked reference/ dir inside an otherwise-stamped install: the
        # doctor reports not-installed, and uninstall leaves everything untouched.
        mktempdir() do dir
            skill = install_agent_guide(dir = dir)
            refdir = joinpath(skill, "reference")
            copy_elsewhere = joinpath(dir, "reference-copy")
            cp(refdir, copy_elsewhere)
            rm(refdir; recursive = true)
            symlink(copy_elsewhere, refdir; dir_target = true)

            @test !agent_guide_status(dir = dir).installed
            @test isempty(uninstall_agent_guide(dir = dir))
            @test isfile(joinpath(skill, "SKILL.md"))
        end

        # (i) Ordering regression: a symlinked .gitignore must be refused BEFORE
        # the skill bundle is written — never a half-installed guide.
        mktempdir() do dir
            external = joinpath(dir, "external-gitignore-order")
            write(external, "external\n")
            symlink(external, joinpath(dir, ".gitignore"))
            @test_throws ArgumentError install_agent_guide(dir = dir)
            @test !ispath(joinpath(dir, ".claude", "skills", "smlma-ecosystem", "SKILL.md"))
        end

        # (j) Same ordering regression for a symlinked AGENTS.md (Codex): GUIDE.md
        # must never appear.
        mktempdir() do dir
            external = joinpath(dir, "external-agents-order")
            write(external, "external\n")
            symlink(external, joinpath(dir, "AGENTS.md"))
            @test_throws ArgumentError install_agent_guide(dir = dir, tool = :codex)
            @test !isfile(joinpath(dir, "smlm-agent-guide", "GUIDE.md"))
        end
    end

    # _expand_home strips repeated leading separators after "~/" too.
    @test SMLMAnalysis._expand_home("~//repo") == joinpath(homedir(), "repo")
    @test SMLMAnalysis._expand_home("~") == homedir()

    # AGENTS.md verbatim restore: content appended after our block (including
    # indentation) must survive an uninstall untouched — not run through strip().
    mktempdir() do dir
        write(joinpath(dir, "AGENTS.md"), "# My project rules\n\nBe careful.\n")
        install_agent_guide(dir = dir, tool = :codex)
        agents_path = joinpath(dir, "AGENTS.md")
        write(agents_path, read(agents_path, String) * "\n    user_code()\n")

        uninstall_agent_guide(dir = dir, tool = :codex)
        final = read(agents_path, String)
        @test occursin("\n    user_code()\n", final)   # indentation intact
        @test occursin("My project rules", final)
        @test !occursin("BEGIN SMLMAnalysis agent-guide", final)
    end

    # `dir = "~/..."` is expanded; it must never create a literal "~" under the cwd,
    # and the doctor's resolved path must actually point under the real home.
    mktempdir() do tmp
        cd(tmp) do
            st = agent_guide_status(dir = "~/smlma-nonexistent-repo-for-test")
            @test !st.installed
            @test st.path == joinpath(
                homedir(), "smlma-nonexistent-repo-for-test",
                ".claude", "skills", "smlma-ecosystem"
            )
            @test !ispath("~")
        end
    end

    # Regression: the OLD (unexpanded) code passed a no-mutation check like the one
    # above too, so assert against a resolved path under a TEMPORARY HOME instead —
    # isolated from the real one — to actually exercise `~` expansion.
    if Sys.isunix()
        mktempdir() do home
            withenv("HOME" => home) do
                mktempdir() do cwd
                    cd(cwd) do
                        p = install_agent_guide(dir = "~/repo")
                        @test startswith(p, joinpath(home, "repo"))
                        @test isfile(joinpath(p, "SKILL.md"))
                        @test !ispath("~")
                    end
                end
            end
        end
    end

    # Claude, track=true → committed (no .gitignore written).
    mktempdir() do dir
        install_agent_guide(dir = dir, track = true)
        @test !isfile(joinpath(dir, ".gitignore"))
    end

    # Codex: stamped bundle + a managed block appended to AGENTS.md that preserves
    # pre-existing content and is idempotent; uninstall strips it back out.
    mktempdir() do dir
        write(joinpath(dir, "AGENTS.md"), "# My project rules\n\nBe careful.\n")
        bundle = install_agent_guide(dir = dir, tool = :codex)
        @test bundle == joinpath(dir, "smlm-agent-guide")
        @test isfile(joinpath(bundle, "GUIDE.md"))
        @test occursin("x-installer: SMLMAnalysis", read(joinpath(bundle, "GUIDE.md"), String))
        @test length(readdir(joinpath(bundle, "reference"))) == 11

        agents = read(joinpath(dir, "AGENTS.md"), String)
        @test occursin("My project rules", agents)                   # user content kept
        @test occursin("BEGIN SMLMAnalysis agent-guide", agents)     # our block added
        @test occursin("smlm-agent-guide/GUIDE.md", agents)

        # Idempotent refresh of our own bundle (no overwrite needed).
        install_agent_guide(dir = dir, tool = :codex)
        agents2 = read(joinpath(dir, "AGENTS.md"), String)
        @test count("BEGIN SMLMAnalysis agent-guide", agents2) == 1  # not duplicated
        @test occursin("My project rules", agents2)

        # Uninstall removes bundle + our block, preserving the user's content.
        removed = uninstall_agent_guide(dir = dir, tool = :codex)
        @test joinpath(dir, "smlm-agent-guide") in removed
        @test !isdir(bundle)
        agents3 = read(joinpath(dir, "AGENTS.md"), String)
        @test occursin("My project rules", agents3)
        @test !occursin("BEGIN SMLMAnalysis agent-guide", agents3)
    end
end

@testset "lab convention conformance" begin
    # Explicit conformance check for the lab skills-installer convention
    # (independent per-package implementations; see the convention doc). Each
    # assert maps to one spec invariant so the lab-guide can cite this block.
    mktempdir() do dir
        skill = install_agent_guide(dir = dir)                    # tool=:claude default
        # 1. Namespaced install dir: <pkgprefix>-<skill>
        @test occursin(r"[/\\]smlma-ecosystem$", skill)
        fm = read(joinpath(skill, "SKILL.md"), String)
        # 2. Provenance stamp: all four x- fields present
        for k in ("x-installer:", "x-source-version:", "x-source-commit:", "x-installed-format:")
            @test occursin(k, fm)
        end
        # 3. copy-never-symlink: installed files are real files, not links
        @test !islink(joinpath(skill, "SKILL.md"))
        # 4. track=false default anchors a .gitignore entry
        @test occursin("/.claude/skills/smlma-ecosystem/", read(joinpath(dir, ".gitignore"), String))
        # 5. own-install refresh is idempotent (no flag, same path)
        @test install_agent_guide(dir = dir) == skill
        # 6. stamp-scoped uninstall removes our own install
        @test uninstall_agent_guide(dir = dir) == [skill]
    end
end

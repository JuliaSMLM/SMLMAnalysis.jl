```@meta
CurrentModule = SMLMAnalysis
```

# AI Assistant Guide

If you work with an AI coding assistant (Claude Code or Codex), SMLMAnalysis can
install a hierarchical, version-stamped guide to the whole JuliaSMLM ecosystem — the
`analyze()` pipeline plus the API reference of every sub-package. The guide is
assembled at install time from the package versions **resolved in your
environment** (each package's `api_overview.md`, or its README as a fallback), so it
never describes an API you do not actually have.

```julia
using SMLMAnalysis
install_agent_guide()                 # Claude Code skill in ./.claude (gitignored by default)
install_agent_guide(track = true)     # …and committed, to share with the repo
install_agent_guide(tool = :codex)    # Codex: AGENTS.md block + reference bundle in this repo
install_agent_guide(scope = :user)    # once for all your projects (~/.claude or ~/.codex)
```

`tool` is `:claude` or `:codex`; `scope` is `:project` (into `dir`, default the
current directory) or `:user` (your home). At project scope the guide is added to
`.gitignore` unless `track = true`.

## What gets installed

| tool | files | discovery |
|---|---|---|
| `:claude` | `.claude/skills/smlma-ecosystem/SKILL.md` + `reference/<Package>.md` per package | Claude Code loads the skill by directory |
| `:codex` | `smlm-agent-guide/GUIDE.md` + `reference/`, plus one managed block in `AGENTS.md` | Codex reads `AGENTS.md`; your existing content is preserved |

Every installed wrapper carries a provenance stamp (`x-installer`,
`x-source-version`, `x-source-commit`, `x-installed-format`) that identifies it as
SMLMAnalysis's own install. Installed copies are never meant to be hand-edited —
re-run the installer to refresh them.

## Refresh, status, uninstall

- **Re-running `install_agent_guide()` refreshes** an install the stamp identifies
  as ours, with no extra flag. A hand-made or foreign directory at the target path
  is **refused** unless you pass `overwrite = true`.
- **`agent_guide_status()`** is the doctor: it reports whether a guide is installed,
  which SMLMAnalysis version and commit produced it, and whether it is `stale`
  relative to the version now resolved.
- **`uninstall_agent_guide()`** removes only what the installer wrote (`SKILL.md` /
  `GUIDE.md`, `reference/`, and the `AGENTS.md` block). The directory itself is
  removed only if nothing else remains, so files you placed alongside are kept.

This installer follows the lab-wide convention for package-shipped assistant guides
(namespaced install directory, `x-` provenance stamp, own-install idempotent refresh,
stamp-scoped uninstall), so it coexists with guides installed by other packages.

## API

```@docs
install_agent_guide
uninstall_agent_guide
agent_guide_status
```

# CellMetPro — Commit Message Guide

Detailed version of the commit conventions summarized in `CONTRIBUTING.md`.

Format: `<type>(<scope>): <short description>`

- Description is lowercase, imperative mood ("add" not "added" / "adds")
- No period at the end
- Max 72 characters on the first line
- Breaking changes: append `!` after type/scope, or add `BREAKING CHANGE:` in footer

---

## Types

| Type       | When to use                                              | Version impact |
|------------|----------------------------------------------------------|----------------|
| `feat`     | A new feature visible to the user or API consumer       | minor bump     |
| `fix`      | A bug fix                                                | patch bump     |
| `perf`     | Performance improvement with no behavior change          | patch bump     |
| `refactor` | Code restructure — no new feature, no bug fix            | none           |
| `style`    | Formatting only — no logic change (whitespace, quotes)   | none           |
| `revert`   | Reverts a previous commit                                | none           |
| `chore`    | Everything that isn't library or app code (see below)    | none           |

A `!` after the type (or `BREAKING CHANGE:` in the footer) triggers a **major bump**:
```
feat(engine)!: change job status values
```

---

## Scopes

| Scope       | What it covers                                                       |
|-------------|----------------------------------------------------------------------|
| `core`      | `cellmetpro/core/` (COMPASS, FBA, preprocessing, progress protocol)  |
| `analysis`  | `cellmetpro/analysis/`                                               |
| `viz`       | `cellmetpro/visualization/`                                          |
| `io`        | `cellmetpro/io.py`, `cellmetpro/data/`                               |
| `report`    | `cellmetpro/reporting/`                                              |
| `cli`       | `cellmetpro/cli.py` (incl. the `app` and `dashboard` subcommands)    |
| `engine`    | `cellmetpro/app/engine/` (projects, files, jobs, runner, tasks, db)  |
| `ui`        | `cellmetpro/app/ui/` (Shiny app and modules)                         |

Scope is optional but strongly encouraged. Omit it only for truly cross-cutting changes.

### `chore` (no scope)

Use plain `chore:` for the launcher, installer (conda recipe, `constructor`, `menuinst`),
tooling (ruff, black, mypy, import-linter, pre-commit), CI workflows, Docker, dependencies,
tests, docs and `context/`, and releases (version bumps, changelog, tags).

---

## Examples

```
feat(core): add progress callback protocol
refactor(core): route compass progress through rich adapter
chore: add app extra with shiny and sqlalchemy
feat(cli): add app subcommand with lazy import
feat(engine): add job submission and startup reconciliation
fix(engine): mark orphaned running jobs as failed on startup
feat(ui): add projects module with trash and restore
chore: reuse running app instance via lock file
chore: test runner with a fast fake task
chore: add windows to ci matrix
feat(cli)!: remove deprecated dashboard subcommand
```

---

## Multi-line commits (when the why needs explaining)

```
fix(engine): mark orphaned running jobs as failed on startup

When the app was quit mid-run, jobs stayed "running" forever because
nothing was left to update them, and the concurrency limit blocked
every new submission. Startup now checks each running job's PID and
fails the ones whose process is gone.
```

The body answers **why**, not what. The diff already shows what changed.

# CellMetPro App — Implementation Roadmap

Single source of truth for the implementation plan of `cellmetpro.app` (see `CONTEXT.md`).
Steps are tackled in order. Check items off as they're completed.

---

## Phase 0 — Library preparation (inside CellMetPro proper)

> The minimum the library needs before an app can sit on top of it.

- [ ] **0.1** — Progress protocol: `cellmetpro/core/progress.py` (callback: step, done, total).
  COMPASS and other long operations report through it. The CLI keeps `rich` through an adapter,
  and `show_progress` behavior stays unchanged for existing users
- [ ] **0.2** — Realistic COMPASS benchmark: genome-scale model, about 5k and 50k cells, record
  runtime and peak memory against `n_processes`. Results set the app's defaults and the
  expectations shown in the UI
- [ ] **0.3** — Packaging: `[app]` extra in `pyproject.toml` (shiny, shinywidgets, sqlalchemy,
  alembic, pydantic), mypy strict override for `cellmetpro.app.*`, import-linter contracts
  (library ↛ app, engine ↛ shiny)
- [ ] **0.4** — `cellmetpro app` CLI subcommand: lazy import, clear message when the extra is
  missing, flags `--host --port --data-dir --open-browser/--no-browser`

---

## Phase 1 — Engine (`cellmetpro/app/engine`, no Shiny)

> The real product. Fully testable without a browser.

- [ ] **1.1** — Config: data dir, host, port, max concurrent jobs. Precedence: CLI flags > env >
  defaults. No hardcoded paths
- [ ] **1.2** — DB layer: sync SQLAlchemy 2.x, SQLite in WAL mode, Alembic, first migration
- [ ] **1.3** — Models ported from `cellmetpro-ui`: Project, File, Job + enums. Job gains
  `params`, `cellmetpro_version`, `pid`, and the `cancelled` status
- [ ] **1.4** — Projects service: CRUD, trash, restore, permanent delete (cascade)
- [ ] **1.5** — Files service: register by path (resolve, exists, extension allowlist),
  `file_type`, status re-check, smart delete, provenance via `job_id`. Never modify inputs
- [ ] **1.6** — Directory listing service for the file browser (rooted, no symlink escape when a
  root is configured)
- [ ] **1.7** — Job lifecycle: submit → queue → subprocess runner → complete/failed/cancelled.
  Concurrency limit, cancel (terminate the process tree), startup reconciliation of orphaned jobs
- [ ] **1.8** — Progress adapter: library progress events → throttled DB writes from the runner
- [ ] **1.9** — Task adapters: parameter model + run function per analysis, starting with COMPASS.
  Outputs go to `<data_dir>/projects/<id>/outputs/`, with a new `.h5ad` written via
  `cellmetpro.io.store_*`
- [ ] **1.10** — Logging: structured, to a log file in the data dir, per-job log capture
- [ ] **1.11** — Engine test suite: services, lifecycle, reconciliation, runner with a fast fake task

---

## Phase 2 — Shiny UI (`cellmetpro/app/ui`)

> A thin client: calls the engine, watches state.

- [ ] **2.1** — Launcher + app skeleton: navbar/sidebar layout, single-instance lock, free-port
  selection, open the browser when the server is ready
- [ ] **2.2** — Projects module: list, create, edit, trash, restore, permanent delete
- [ ] **2.3** — File browser + registration module (server-side listing, type selection, status
  badges)
- [ ] **2.4** — Jobs panel: `reactive.poll` on job state, progress bars, cancel, view log, notice
  that closing the tab doesn't stop jobs
- [ ] **2.5** — COMPASS form → submit job (validated parameters, sensible defaults from 0.2)
- [ ] **2.6** — Results viewer: Plotly via `shinywidgets`, reusing
  `cellmetpro.visualization.interactive`. Check responsiveness with 50k+ cells (WebGL traces)
- [ ] **2.7** — Remaining analyses, one at a time (form + task + view): differential (incl.
  pseudo-bulk), clustering, pathway, trajectory, perturbation, batch correction
- [ ] **2.8** — Report export (`cellmetpro.reporting`)
- [ ] **2.9** — Settings page: data dir, `n_processes`, solver, concurrency
- [ ] **2.10** — First-run experience: bundled sample dataset (`cellmetpro.data`), guided demo
  project
- [ ] **2.11** — UI tests (`shiny.playwright` + `pytest-playwright`) and accessibility pass
- [ ] **2.12** — Deprecate the Streamlit dashboard (warning in `cellmetpro dashboard`, removal
  planned for the next minor release)

---

## Phase 3 — Distribution

> From a `pip` extra to a website download button.

- [ ] **3.1** — CI matrix (Linux, macOS, Windows) including the `[app]` extra and Playwright tests
- [ ] **3.2** — Conda package for CellMetPro (`noarch: python`), built in CI into a local channel.
  Verify conda-forge `cobra`/`swiglpk` on osx-arm64
- [ ] **3.3** — `constructor` installers for macOS (`.pkg`), Windows (`.exe`) and Linux (`.sh`),
  with `menuinst` shortcuts running `cellmetpro app --open-browser`, built on tag
- [ ] **3.4** — Signing: macOS notarization. Decide on Windows code signing, or document the
  SmartScreen warning
- [ ] **3.5** — Download website (GitHub Pages): OS-detected download button linking to GitHub
  Releases assets, quick start, screenshots
- [ ] **3.6** — Release workflow: tag → PyPI + installers + GitHub Release + changelog entry
- [ ] **3.7** — Submit CellMetPro to conda-forge (staged-recipes)

---

## Phase 4 — Hosted deployment (optional, later)

- [ ] **4.1** — Docker image running `cellmetpro app --host 0.0.0.0 --no-browser --data-dir /data`
- [ ] **4.2** — Hardening for non-localhost use: file browser confined to the data root, auth via
  reverse proxy or ShinyProxy (one container per user)
- [ ] **4.3** — Public demo instance on the sample dataset

---

## Phase 5 — Community

- [ ] **5.1** — Zenodo DOI per GitHub release
- [ ] **5.2** — bio.tools registry entry
- [ ] **5.3** — Docs site + short demo video
- [ ] **5.4** — JOSS paper once stable

---

## Legend

| Symbol | Meaning |
|--------|---------|
| `[ ]`  | Not started |
| `[~]`  | In progress |
| `[x]`  | Done |

# CellMetPro App — Project Context

> Supersedes the original CellMetPro UI plan (FastAPI + React + Electron, separate repo).
> That version is preserved in the git history of `omicspilot/cellmetpro-ui`.
> Decision date: 2026-10-09. Moved into the CellMetPro repo (`context/`) the same day; work
> happens on the `feature-add-shiny-app` branch.

## What this is

A graphical app for **CellMetPro** (https://github.com/omicspilot/CellMetPro), the Python package
for single-cell RNA-seq metabolic analysis (COMPASS scoring, FBA, differential, clustering, pathway,
trajectory, perturbation, reporting). CellMetPro is on PyPI with a CLI and a Python API.

**The goal is to make CellMetPro usable by scientists with little programming or command-line
experience.** They download an installer from a website, double-click a "CellMetPro" shortcut, and
the app opens in their web browser. They never see a terminal, `pip`, or Python.

The app lives **inside the CellMetPro package** as an optional subpackage, `cellmetpro.app`,
installed with `pip install cellmetpro[app]` and launched with `cellmetpro app`. It replaces the
existing Streamlit dashboard.

Production quality, open source (MIT), cross-platform (macOS, Windows, Linux).

---

## Decisions and why

These were evaluated on 2026-10-09. Don't reopen them without new information.

| Decision | Why |
|---|---|
| **The primary goal is shipping to non-CLI scientists**, ahead of learning a particular tech stack | Stated by the developer. Every choice below follows from it |
| **No Electron, no desktop shell.** The browser is the window | The "app feel" isn't important. What matters is downloading from a website. Electron's window, IPC, process manager and auto-updater add a large surface area for little value to this audience |
| **There is still an installer** | Something has to start the Python process, and it can only be the user in a terminal (ruled out), the user via a double-click, or someone else (hosting). A local app without a CLI therefore needs an installer that bundles Python and all dependencies |
| **The installer is built with conda `constructor`**, with `menuinst` shortcuts | Standard in scientific Python (napari ships its installers this way). It's a config file plus a CI job, not a sub-project |
| **Shiny for Python** for the UI | One language, no frontend build step, direct library calls. With Electron gone, React's remaining advantage (UI polish) doesn't justify a second stack for a solo developer. The reactive model carries over to R Shiny, which is widespread in bioinformatics |
| **Merged into the CellMetPro package**, not a separate repo or package | The old plan already pinned server version = library version, so they were one release unit anyway. Calling the library directly removes the HTTP, schema and generated-client layers. Merging costs nothing at import time (see "Package weight") |
| **No remote mode** (desktop client talking to a remote server) | That client/server split was the main source of complexity. A *hosted* deployment (the same app running on a lab server, opened by URL) stays possible later at almost no cost, because it's the same artifact |
| **Single user per running instance** | The engine has no notion of identity. Multi-user hosting, if ever needed, means one container per user (ShinyProxy model), not multi-tenancy in the data model |
| **Streamlit dashboard is deprecated** | One UI only. `cellmetpro dashboard` keeps working for one minor release with a deprecation warning pointing to `cellmetpro app`, then it's removed along with the `dashboard` extra |

### Package weight

Merging doesn't make `import cellmetpro` heavier. Measured on 2026-10-09: `import cellmetpro`
takes about 0.1 s and 19 MB of memory, and loads none of scanpy, plotly, matplotlib, cobra or
sklearn, because `cellmetpro/__init__.py` loads subpackages lazily through `__getattr__`. Shiny and
its web dependencies install only with the `[app]` extra.

**Rules that keep it this way:**
- Nothing outside `cellmetpro/app/` imports `cellmetpro.app`.
- `app` is not listed in the lazy `__all__` in `cellmetpro/__init__.py`.
- The `cellmetpro app` CLI handler imports the app *inside* the function, as `run_dashboard`
  already does for Streamlit, and prints a clear "install `cellmetpro[app]`" message on `ImportError`.

---

## Architecture

### Layering

```
cellmetpro/core, analysis, visualization, io, ...   ← the library (unchanged role)
        ▲
cellmetpro/app/engine/    pure Python, NO Shiny imports   ← projects, files, jobs, runner
        ▲
cellmetpro/app/ui/        Shiny                          ← replaceable client
```

- Dependencies point one way: `ui → engine → library`. The library never imports `app`, and the
  engine never imports Shiny.
- The UI never runs heavy computation or spawns processes itself. It calls the engine and watches
  state.
- Enforce the layering mechanically, for example with `import-linter` contracts in CI, not by
  convention alone.

### Proposed layout

```
cellmetpro/
├── core/ analysis/ visualization/ reporting/ data/ models/ io.py cli.py   # existing
└── app/
    ├── __init__.py          # imports nothing heavy
    ├── launcher.py          # serve(host, port, data_dir, open_browser): single-instance lock,
    │                        # free-port selection, start server, open browser when ready
    ├── engine/
    │   ├── config.py        # data dir, host, port, max concurrent jobs (CLI flags > env > defaults)
    │   ├── db.py            # SQLAlchemy engine/session, SQLite in WAL mode
    │   ├── models.py        # Project, File, Job + enums (ported from cellmetpro-ui)
    │   ├── migrations/      # Alembic
    │   ├── projects.py      # service: CRUD, trash, restore, permanent delete
    │   ├── files.py         # service: register by path, type, status checks, provenance
    │   ├── jobs.py          # service: submit, queue, cancel, list, reconcile on startup
    │   ├── runner.py        # entry point of the job subprocess: load job → run task → write results
    │   ├── progress.py      # adapter: library progress events → throttled DB writes
    │   └── tasks/           # one adapter per analysis: params model + run function
    │       ├── compass.py
    │       └── ...
    └── ui/
        ├── app.py           # Shiny App assembly
        └── modules/         # one Shiny module per feature (projects, files, jobs, compass, ...)
```

### Long-running jobs

COMPASS solves one linear program per cell with GLPK, spread over a `ProcessPoolExecutor`.
Real runs on genome-scale models can take hours. Therefore:

- **Jobs never run inside a Shiny session.** Shiny's `ExtendedTask` is tied to a browser session,
  which ends when the tab closes.
- **Each job runs in its own subprocess** (`python -m cellmetpro.app.engine.runner <job_id>`) started
  by the engine. It must not be a daemonic `multiprocessing.Process`, because daemonic processes
  can't create the child processes COMPASS's process pool needs.
- **SQLite is the single source of truth** for job state (`pending → running → complete | failed |
  cancelled`), progress, PID, and timestamps. WAL mode lets the runner write while the UI reads.
  Progress writes are throttled (about once per second at most).
- **The UI watches with `reactive.poll`**, using a cheap check such as the latest `updated_at` of
  the project's jobs.
- **Concurrency limit:** one running job at a time by default, configurable, because COMPASS
  already uses all available cores.
- **Startup reconciliation:** when the app starts, any job marked `running` whose PID is gone is
  marked `failed` with a clear message.
- **Closing the browser tab doesn't stop jobs.** Quitting the app (stopping the server process)
  does. The UI says so.

### Progress reporting (library change)

Today the library's progress output is hard-wired to `rich` terminal bars
(`core/compass.py`, `show_progress`). Add a small **library-level progress protocol**
(for example `cellmetpro/core/progress.py`: a callback receiving step name, done, and total). The
CLI keeps `rich` through an adapter, and the app writes to the DB through another. This is the only
change the app requires inside the library proper.

### Files

- **Register by path, never upload.** The app runs on the user's machine, so files stay where they
  are. The DB stores the absolute path and metadata.
- **File browser on the server side.** The app lists directories itself instead of using the
  browser's file input, which would copy multi-GB `.h5ad` files. The same browser works for a later
  hosted deployment, rooted at the mounted data directory.
- Accepted inputs: `.h5ad`, `.csv`, `.mtx` (plus pre-computed result tables).
- `file_type` enum and `status` (`available`, `missing`, ...) as in the old design. Status is
  re-checked when a project is opened.
- **Never modify the user's input files.** Outputs go to `<data_dir>/projects/<project_id>/outputs/`.
  When results are written into AnnData (`cellmetpro.io.store_*`), they go into a *new* `.h5ad`
  in the outputs folder.
- Output files carry `job_id`, so every result traces back to job → parameters → input files →
  CellMetPro version.

### Data model (ported from cellmetpro-ui)

- **Project:** name, description, optional `workspace_path`, timestamps, `deleted_at` for trash.
  Trash, restore, and permanent delete (cascades to files and jobs, never deletes registered
  originals).
- **File:** project, filename, size, path, `file_type`, `status`, optional `job_id`.
- **Job:** project, `analysis_type`, status, step, progress, message, parameters (JSON), results
  metadata (JSON), CellMetPro version, PID, timestamps.
- Analyses are modular and independent. No forced pipeline order; a user can register a
  pre-computed differential table and visualize it directly.

### Launcher behavior

- Binds to `127.0.0.1` by default. Binding anywhere else requires an explicit `--host` and logs a
  warning, because the app has no auth and the file browser exposes the filesystem.
- Picks a free port if the default is taken.
- **Single instance:** a lock file in the data directory records the running port. A second launch
  (double-clicking the shortcut again) opens the existing instance instead of starting a second
  server against the same SQLite DB.
- Opens the default browser only once the server responds.
- Default data dir is `~/.cellmetpro/`; everything can be overridden (`--data-dir`, env). No
  hardcoded paths.

---

## Starting point in this repo

What the plan touches in CellMetPro as it stands today (verified 2026-10-09, v0.2.0).
Line numbers drift; re-check before relying on them.

- **Progress (0.1):** `rich.progress` is imported in `cellmetpro/core/compass.py` and used behind
  `self.config.show_progress` in four places (COMPASS scoring loops). `rich` is also used by
  `cli.py` (via `rich-argparse`). No other library module reports progress yet.
- **Lazy loading:** `cellmetpro/__init__.py` exposes `core, analysis, visualization, data,
  reporting` through `__getattr__`. Keep `app` out of that list.
- **CLI pattern to copy (0.4):** `run_dashboard` in `cellmetpro/cli.py` already imports Streamlit
  inside the handler and prints an install hint on `ImportError`. `cellmetpro app` follows the
  same shape; subcommands are dispatched in `main()`.
- **Streamlit surface to deprecate (2.12):** `dashboard` subcommand + `run_dashboard` in
  `cli.py`, `cellmetpro/visualization/dashboard.py`, the `dashboard` extra (also referenced by
  `all`), `streamlit` in `conda_config.yaml`, the coverage `omit` entry for `dashboard.py`, and
  README mentions.
- **Packaging (0.3):** `pyproject.toml` has extras `dev`, `dashboard`, `seurat`, `all`. Add `app`,
  and decide whether `all` includes it. mypy is global and lenient
  (`ignore_missing_imports = true`, not strict), so strictness for `cellmetpro.app.*` goes in a
  `[[tool.mypy.overrides]]` block.
- **Package discovery:** `[tool.setuptools.packages.find]` includes `cellmetpro*`, so
  `cellmetpro.app` subpackages are picked up automatically. Non-Python files (Alembic
  `script.py.mako`, `alembic.ini`, static assets) need `[tool.setuptools.package-data]` (or a
  `MANIFEST.in`) or they won't ship in the wheel.
- **`.gitignore` traps:** it ignores `data/`, `results/`, `output/`, `**/cache/`, `*.log`,
  `*.html` and `.claude`. Avoid those names for directories under `cellmetpro/app/` (or add
  explicit `!` exceptions), and keep test fixtures out of them.
- **CI (3.1):** `.github/workflows/ci.yml` runs on push/PR to `main` and `dev`, matrix
  `ubuntu-latest` + `macos-latest` × Python 3.10–3.12, GLPK installed per OS. No Windows yet, and
  no `[app]` extra or Playwright job.
- **Conda (3.2):** `conda_config.yaml` is a dev environment, not a recipe. It pip-installs
  `cobra>=0.29` and `pyarrow` and includes `streamlit`; the installer needs a real conda recipe.
- **Branch flow:** feature branch → `dev` → `main` (see `CONTRIBUTING.md`). Commit messages
  follow `context/COMMIT_GUIDE.md`.
- **Checks before every commit:** black, ruff, mypy, pytest (as listed in `CONTRIBUTING.md`).

---

## Distribution

| Channel | Audience | How |
|---|---|---|
| **Installer** (primary) | Scientists without CLI experience | Download page → `.pkg` (macOS), `.exe` (Windows), `.sh` (Linux) built with conda `constructor`. Bundles Python + `cellmetpro[app]`. A `menuinst` shortcut runs `cellmetpro app --open-browser` |
| **PyPI** | Developers, notebook users | `pip install cellmetpro` (library + CLI) or `cellmetpro[app]` |
| **Docker** (later) | Labs, core facilities | Same `cellmetpro app` command, `--host 0.0.0.0 --data-dir /data`, one container per user, auth via reverse proxy |

Facts to plan around:
- `constructor` bundles **conda packages only**. CellMetPro and all its dependencies must be
  available as conda packages. CellMetPro is pure Python (`noarch: python`). **Bioconda doesn't
  build for Windows**, so the long-term home is **conda-forge**. Short term, CI can build a local
  noarch package and feed it to `constructor`.
- `conda_config.yaml` currently installs `cobra` via pip for Apple Silicon. Verify that the
  conda-forge `cobra` and `swiglpk` work on osx-arm64 before relying on them in the installer.
- Unsigned installers trigger macOS Gatekeeper and Windows SmartScreen warnings. macOS
  notarization requires an Apple Developer account. Windows signing can come later, with the
  warning documented on the download page.
- The installer will be a few hundred MB (Python + scientific stack). That's expected.
- Reference models and large data are never bundled beyond the small sample dataset already in
  `cellmetpro.data`.
- **Shinylive is not an option:** cobra, GLPK and process pools can't run in WebAssembly.

Community visibility (once stable): Zenodo DOI per release, bio.tools entry, JOSS paper, docs
site with a short demo video, and a first-run tutorial on the bundled sample data.

---

## Quality standards

Same bar as the original plan, adapted to the stack:
- **mypy `--strict` for `cellmetpro.app.*`** (via a per-module override; the rest of CellMetPro
  keeps its current settings).
- **ruff + black**, consistent with CellMetPro.
- **The engine is tested without Shiny** (pytest): services, job lifecycle, reconciliation, runner
  with a fast fake task.
- **UI tests** with Shiny's Playwright controllers (`shiny.playwright`) and `pytest-playwright`.
- **Layering enforced in CI** (import-linter or equivalent).
- **CI matrix** on Linux, macOS, Windows, including the `[app]` extra and an installer smoke test
  on release.
- No `print` in production code; structured logging. The log file lives in the data dir and the
  UI can show it.
- Keep a Changelog + semantic versioning (one version for library, CLI and app).
- Accessibility: keyboard navigation, labels on every input, ARIA on custom components.

---

## Porting from `cellmetpro-ui`

Source paths below are relative to the `cellmetpro-ui` repo (checked out next to this one, at
`../cellmetpro-ui`).

**Bring over (adapt to sync SQLAlchemy; there's no reason for async at single-user scale):**
- `packages/server/cellmetpro_server/models.py`: enums (`FileType`, `FileStatus`, `AnalysisType`),
  Project, File, Job.
- Project trash, restore and permanent-delete semantics (`routers/projects.py`).
- File registration logic (`routers/io.py`): path resolution, existence check, extension
  allowlist, smart delete (never delete registered originals).
- Error codes (`utils/const.py`) → engine exceptions.
- Config defaults (`config.py`): data dir under `~/.cellmetpro/`.

**Leave behind:** FastAPI routers, WebSocket progress, multipart upload, OpenAPI export, the
whole React/Vite/Zustand/TanStack/openapi-ts frontend, the Electron plan.

The `cellmetpro-ui` repo is archived once the port is done.

---

## The developer and how to collaborate

**Oumar Ndiaye**, bioinformatics engineer and author of CellMetPro. Advanced Python (built
CellMetPro end-to-end with CI/CD, 80% coverage, PyPI). Strong JavaScript (React, Vue). Not a fan
of Java.

Learning through this project: Shiny's reactive model, process and job management, SQLite
concurrency, conda packaging and `constructor` installers, release pipelines, later Docker.

Collaboration preferences:
- **No code unless explicitly asked** ("give me the code", "write it", "implement it").
- **Theory first.** For each step: an overview that weaves need and solution together, a list of
  references (official docs, quality deep dives), then progressive pseudocode skeletons with
  inline comments explaining the *why*: what breaks if a line is skipped, which alternative was
  rejected and why. The developer writes the real code.
- **Hint-based corrections**, not rewrites.
- **Scope discipline:** change only what's explicitly requested. No proactive cleanup of adjacent
  files.
- Production quality, no "clean this up later".

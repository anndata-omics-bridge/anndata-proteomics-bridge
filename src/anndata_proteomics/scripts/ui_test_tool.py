#!/usr/bin/env python3
"""APB test-data browser → convert (background job) → inspect (marimo GUI).

Browse the ProteoBench test corpus, filter by software / size / conversion target, and convert a
selected dataset to an AnnData level or a MuData. The conversion runs as a **background
subprocess** (anndata_proteomics.scripts.convert_one) streaming to a log file; the UI polls with
mo.ui.refresh and shows a live status panel — so heavy fragment/MuData runs don't freeze the app
and the log/result open in a browser. See the marimo-background-jobs skill + TODO/HOWTO/test_gui.md.

Run via the top-level Makefile (`make ui`) or directly:
    marimo run src/anndata_proteomics/scripts/ui_test_tool.py
"""

import marimo

app = marimo.App(width="full")


@app.cell
def _():
    import sys
    from datetime import datetime

    import marimo as mo

    from anndata_proteomics.scripts import _ui_panels as panels
    from anndata_proteomics.scripts import _ui_support as ui
    from anndata_proteomics.scripts import jobrunner as runner

    return datetime, mo, panels, runner, sys, ui


@app.cell
def _(mo):
    mo.md(
        """
        # APB test-data browser

        Pick a **target**, filter the catalog, select a dataset, and **Convert**. The job runs in
        the background; the panel on the right shows live status + console log, with the result and
        log openable in the browser.
        """
    )
    return


@app.cell
def _(ui):
    catalog = ui.load_catalog()
    return (catalog,)


@app.cell
def _(catalog, mo, ui):
    targets = ui.LEVELS + [ui.MUDATA]
    softwares = ["All"] + sorted(catalog["software_name"].unique()) if len(catalog) else ["All"]
    max_mb = float(catalog["size_mb"].max()) if len(catalog) else 1000.0

    target_dd = mo.ui.dropdown(options=targets, value="ion", label="Target")
    software_dd = mo.ui.dropdown(options=softwares, value="All", label="Software")
    size_slider = mo.ui.slider(
        start=0, stop=max_mb, value=max_mb, step=1, label="Max size (MB)", show_value=True
    )
    return size_slider, software_dd, target_dd


@app.cell
def _(mo):
    converted_refresh = mo.ui.run_button(label="Refresh outputs")
    return (converted_refresh,)


@app.cell
def _(mo):
    # Background-job state (see the marimo-background-jobs skill):
    #   get_job   – the jobrunner.Job handle for the current/last conversion
    #   get_server– the static HTTP server serving a finished job's output folder
    #   get_done  – run_key of a job already observed finished (stops auto-poll cleanly)
    get_job, set_job = mo.state(None)
    get_server, set_server = mo.state(None)
    get_done, set_done = mo.state(None)
    convert_button = mo.ui.run_button(label="Convert ▶", kind="success")
    return convert_button, get_done, get_job, get_server, set_done, set_job, set_server


@app.cell
def _(get_done, get_job, mo):
    # Auto-poll only while a job runs. Read the state getters (not graph vars) so this cell stays
    # outside the run cell's dependency cycle; depending on get_done() rebuilds it (auto-poll off)
    # once the process exits.
    get_done()
    _job = get_job()
    _running = _job is not None and _job.process.poll() is None
    log_refresh = mo.ui.refresh(
        options=["1s", "2s", "5s", "10s"],
        default_interval="2s" if _running else None,
        label="Refresh log",
    )
    return (log_refresh,)


@app.cell
def _(catalog, size_slider, software_dd, target_dd, ui):
    filtered = ui.filter_catalog(
        catalog,
        target=target_dd.value,
        software=software_dd.value,
        max_size_mb=size_slider.value,
    )
    return (filtered,)


@app.cell
def _(filtered, mo):
    # Include targets_str + param_path so the selected row carries everything the Convert needs.
    _cols = [
        "software_name",
        "software_version",
        "nr_prec",
        "size_mb",
        "slug",
        "targets_str",
        "param_path",
        "input_file_path",
    ]
    left_table = mo.ui.table(
        filtered[_cols] if len(filtered) else filtered,
        selection="single",
        page_size=6,
        label=f"Datasets ({len(filtered)})",
    )
    return (left_table,)


@app.cell
def _(
    convert_button,
    datetime,
    get_done,
    get_job,
    left_table,
    log_refresh,
    mo,
    panels,
    runner,
    set_done,
    set_job,
    set_server,
    sys,
    target_dd,
    ui,
    get_server,
):
    job = get_job()
    running = job is not None and job.process.poll() is None

    _sel = left_table.value
    _row = _sel.iloc[0] if (_sel is not None and len(_sel)) else None
    _target = target_dd.value

    # Start a new background conversion on click (replacing any previous job). Each run writes to
    # its own output dir, so distinct conversions never clobber each other.
    if convert_button.value and _row is not None and not running:
        if job is not None:
            runner.terminate_job(job)
        _run_key = runner.make_run_key(_row["input_file_path"], _row["slug"], _target)
        _stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
        _outdir = ui.CONVERTED_DIR / f"{_stamp}_{_row['slug']}_{_target}"
        _cmd = [
            sys.executable,
            "-m",
            "anndata_proteomics.scripts.convert_one",
            "--input",
            _row["input_file_path"],
            "--slug",
            _row["slug"],
            "--target",
            _target,
            "--params",
            _row["param_path"],
            "--outdir",
            str(_outdir),
        ]
        job = runner.start_job(_cmd, _outdir, log_file=_outdir / "console.log", run_key=_run_key)
        set_job(job)
        running = True

    if job is None:
        status_panel = mo.md("*Select a dataset + target, then click **Convert**.*")
    else:
        _ = log_refresh.value  # depend on the refresh tick → re-inspect each poll
        status = runner.inspect_job(job, report_glob="result.*")
        # Flag the finished job once (run_key-guarded) so the refresh widget rebuilds with
        # auto-poll off without looping.
        if not status.running and get_done() != job.run_key:
            set_done(job.run_key)
        _cmd = list(job.command)
        _job_slug = _cmd[_cmd.index("--slug") + 1] if "--slug" in _cmd else ""
        _job_target = _cmd[_cmd.index("--target") + 1] if "--target" in _cmd else ""
        status_panel = panels.build_status_panel(
            status,
            refresh_widget=log_refresh,
            get_server=get_server,
            set_server=set_server,
            dataset_label=_job_slug,
            target=_job_target,
        )
    return (status_panel,)


@app.cell
def _(converted_refresh, get_done, ui):
    get_done()
    _ = converted_refresh.value
    converted_runs = ui.list_converted_runs()
    return (converted_runs,)


@app.cell
def _(converted_runs, mo, ui):
    converted_display = ui.converted_runs_table(converted_runs)
    converted_table = mo.ui.table(
        converted_display,
        selection="single",
        page_size=5,
        label=f"Converted outputs ({len(converted_display)})",
    )
    return (converted_table,)


@app.cell
def _(converted_runs, converted_table, mo, panels, ui):
    _sel = converted_table.value
    _selected = _sel.iloc[0] if (_sel is not None and len(_sel)) else None
    _row = None
    if _selected is not None and not converted_runs.empty:
        _matches = converted_runs[converted_runs["run_name"] == _selected["run_name"]]
        _row = _matches.iloc[0] if len(_matches) else None
    if _row is None:
        result_viewer = mo.md("*Select a converted output to inspect it.*")
    elif not _row["result_path"]:
        result_viewer = mo.md(f"**No result file for selected run.** Status: `{_row['status']}`")
    else:
        try:
            _obj = ui.load_converted_result(_row["result_path"])
            _summary = ui.summarize(_obj)
            result_viewer = mo.vstack(
                [
                    mo.md(f"**{_row['result_type']} result**  \n`{_row['result_path']}`"),
                    panels.build_summary_panel(_summary),
                ]
            )
        except Exception as exc:  # noqa: BLE001 - render load failures in the GUI.
            result_viewer = mo.md(f"**Could not load result:** `{exc}`")
    return (result_viewer,)


@app.cell
def _(
    convert_button,
    converted_refresh,
    converted_table,
    left_table,
    mo,
    result_viewer,
    size_slider,
    software_dd,
    status_panel,
    target_dd,
):
    input_panel = mo.vstack(
        [
            mo.hstack([target_dd, software_dd, size_slider], justify="start", gap=2),
            left_table,
            convert_button,
            status_panel,
        ]
    )
    converted_panel = mo.vstack(
        [
            mo.hstack([mo.md("## Converted outputs"), converted_refresh], justify="space-between"),
            converted_table,
            mo.md("## Result viewer"),
            result_viewer,
        ]
    )
    mo.vstack([input_panel, converted_panel], gap=2)
    return


if __name__ == "__main__":
    app.run()

"""marimo status-panel for background conversion jobs (the only marimo-aware UI glue).

Renders a Running / Finished / Failed callout with the live console-log tail, the log-file path,
a results download, and an "Open output folder" link served over a real static HTTP server (so the
log and result open in a browser tab — not a file:// or data: link). Mirrors MetaboStatsHub's
apps/job_status.build_run_status_panel; see the marimo-background-jobs skill.
"""

from __future__ import annotations

from collections.abc import Callable

import marimo as mo

from anndata_proteomics.scripts import jobrunner as runner


def build_status_panel(
    status: runner.JobStatus,
    *,
    refresh_widget: object,
    get_server: Callable[[], runner.StaticServer | None],
    set_server: Callable[[runner.StaticServer | None], None],
    dataset_label: str = "",
    target: str = "",
    max_log_chars: int = 6000,
) -> mo.Html:
    """Build the job-status callout. ``status`` comes from ``jobrunner.inspect_job``."""
    if status.running:
        state_label, kind = "Running…", "info"
    elif status.success:
        state_label, kind = "Finished (exit 0)", "success"
    else:
        state_label, kind = f"Failed (exit {status.returncode})", "danger"

    header = f"{dataset_label} → {target}" if dataset_label else "Conversion"
    items: list[object] = [
        mo.md(f"**{header}**"),
        mo.md(f"**Status:** {state_label}"),
        mo.md(f"Output: `{status.output_dir}`"),
        mo.md(f"Log: `{status.log_file}`"),
    ]

    if not status.running and status.output_dir.exists():
        # Serve the run's folder over HTTP so the console.log and result.* open in a browser
        # tab. The .zip download stays available even on failure so the log is always grabbable.
        server = get_server()
        if (
            server is None
            or not server.running
            or server.root_dir != status.output_dir.resolve()
        ):
            if server is not None and server.running:
                server.process.terminate()
            server = runner.start_static_server(status.output_dir)
            set_server(server)
        items.append(mo.md(f"📁 **[Open output folder ↗]({server.base_url}/)**"))
        items.append(
            mo.download(
                data=lambda: runner.zip_dir_to_bytes(status.output_dir),
                filename=f"{(dataset_label or 'result').replace(' ', '_')}_{target}.zip",
                label="Download (.zip)",
            )
        )

    if status.running:
        items.append(refresh_widget)

    log_text = status.log_text or "*(no log output yet)*"
    items.append(
        mo.accordion({"Console log": mo.md(f"```\n{log_text[-max_log_chars:]}\n```")}, lazy=False)
    )
    return mo.callout(mo.vstack(items, gap=0.4), kind=kind)

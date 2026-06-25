"""marimo status-panel for background conversion jobs (the only marimo-aware UI glue).

Renders a Running / Finished / Failed callout with the live console-log tail, the log-file path,
a results download, and an "Open output folder" link served over a real static HTTP server (so the
log and result open in a browser tab — not a file:// or data: link). Mirrors MetaboStatsHub's
apps/job_status.build_run_status_panel; see the marimo-background-jobs skill.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

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
        if server is None or not server.running or server.root_dir != status.output_dir.resolve():
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


# --- Result-summary panel ----------------------------------------------------
# Renders the dict from ``_ui_support.summarize`` (MuData or AnnData) as a
# collapsible, tabular layout instead of a raw JSON code block. Stat badges up
# top, search parameters and each modality as accordion sections, and the raw
# JSON kept as a final section via marimo's built-in tree viewer (``mo.json``).

_MOD_KEYS = ("fixed_mods", "variable_mods")


def _fmt_value(value: Any) -> str:
    """Flatten a search-parameter value to a readable one-liner for a table cell."""
    if isinstance(value, bool):
        return "✓ yes" if value else "✗ no"
    if isinstance(value, Mapping):
        if set(value) == {"value"}:
            return _fmt_value(value["value"])
        if "label" in value:
            mode = value.get("mode")
            return f"{value['label']}" + (f" ({mode})" if mode else "")
        return ", ".join(f"{k}: {_fmt_value(v)}" for k, v in value.items())
    if isinstance(value, (list, tuple)):
        return ", ".join(_fmt_value(v) for v in value) if value else "—"
    return str(value)


def _kv_table_md(data: Mapping[str, Any]) -> str:
    """A two-column key/value markdown table (keys in bold)."""
    if not data:
        return ""
    rows = "\n".join(f"| **{key}** | {_fmt_value(val)} |" for key, val in data.items())
    return f"| parameter | value |\n|---|---|\n{rows}"


def _mods_table_md(name: str, mods: Sequence[Mapping[str, Any]]) -> str:
    """A bold-titled table for a modification list (fixed_mods / variable_mods)."""
    title = f"**{name}** ({len(mods)})"
    if not mods:
        return f"{title} — *none*"
    rows = "\n".join(
        f"| `{m.get('name', '')}` | {m.get('position', '')} | {m.get('mod_type', '')} |"
        for m in mods
    )
    return f"{title}\n\n| modification | position | type |\n|---|---|---|\n{rows}"


def _list_md(name: str, items: Sequence[Any]) -> str:
    """A bold list name followed by its items as wrapping inline-code chips."""
    title = f"**{name}** ({len(items)})"
    if not items:
        return f"{title} — *none*"
    chips = " ".join(f"`{item}`" for item in items)
    return f"{title}<br>{chips}"


def _params_panel(params: Mapping[str, Any]) -> mo.Html:
    """Search parameters: scalar key/value table + modification tables + leftover lists."""
    scalars = {k: v for k, v in params.items() if k not in _MOD_KEYS and not isinstance(v, list)}
    list_keys = [k for k, v in params.items() if k not in _MOD_KEYS and isinstance(v, list)]
    blocks = [_kv_table_md(scalars)]
    blocks += [_mods_table_md(k, params[k]) for k in _MOD_KEYS if k in params]
    blocks += [_list_md(k, params[k]) for k in list_keys]
    return mo.md("\n\n".join(block for block in blocks if block))


def _software_badge(summary: Mapping[str, Any]) -> mo.Html | None:
    """A `software` badge (name + version caption) from search parameters, or None."""
    params = summary.get("search_parameters") or {}
    name = params.get("software_name")
    if not name:
        return None
    return mo.stat(name, label="software", caption=params.get("software_version"), bordered=True)


def _shape_row(summary: Mapping[str, Any], *, lead: list[object] | None = None) -> mo.Html:
    """kind · observations · features badges for an AnnData / modality (optional ``lead`` badges)."""
    n_obs, n_vars = summary.get("shape", ("?", "?"))
    feat = f"{n_vars:,}" if isinstance(n_vars, int) else n_vars
    badges = list(lead or [])
    badges += [
        mo.stat(summary.get("kind", "AnnData"), label="kind", bordered=True),
        mo.stat(n_obs, label="observations", bordered=True),
        mo.stat(feat, label="features", bordered=True),
    ]
    return mo.hstack(badges, justify="start", gap=1, wrap=True)


def _anndata_panel(
    summary: Mapping[str, Any], *, include_params: bool = True, show_shape: bool = True
) -> mo.Html:
    """Body for one AnnData (or one MuData modality): shape badges then column/layer lists."""
    items: list[object] = []
    if show_shape:
        items.append(_shape_row(summary))
    list_md = "\n\n".join(
        _list_md(name, summary[name])
        for name in ("obs_columns", "var_columns", "layers", "uns_keys")
        if name in summary
    )
    if list_md:
        items.append(mo.md(list_md))
    if include_params and summary.get("search_parameters"):
        items.append(_params_panel(summary["search_parameters"]))
    return mo.vstack(items, gap=0.6)


def build_summary_panel(summary: Mapping[str, Any]) -> mo.Html:
    """Render a ``_ui_support.summarize`` dict (MuData or AnnData) as a collapsible panel.

    Stat badges up top; search parameters and each modality as accordion sections; the
    raw dict kept as a final ``mo.json`` tree section so nothing is hidden.
    """
    if summary.get("kind") == "MuData":
        modalities: dict[str, Any] = summary.get("modalities", {})
        software = _software_badge(summary)
        header_items: list[object] = [software] if software else []
        header_items += [
            mo.stat("MuData", label="kind", bordered=True),
            mo.stat(summary.get("n_obs", "—"), label="observations", bordered=True),
            mo.stat(len(modalities), label="modalities", bordered=True),
        ]
        # One feature count per modality — a summed "total" mixes incomparable axes
        # (ion precursors vs protein groups), so report each modality separately.
        for name, ad in modalities.items():
            n_vars = ad.get("shape", ("?", "?"))[1]
            feat = f"{n_vars:,}" if isinstance(n_vars, int) else n_vars
            header_items.append(mo.stat(feat, label=name, caption="features", bordered=True))
        header = mo.hstack(header_items, justify="start", gap=1, wrap=True)
        sections: dict[str, object] = {}
        if summary.get("search_parameters"):
            sections["🔬 Search parameters"] = _params_panel(summary["search_parameters"])
        for name, ad in modalities.items():
            n_obs, n_vars = ad.get("shape", ("?", "?"))
            sections[f"🧬 {name} — {n_obs} × {n_vars}"] = _anndata_panel(ad)
        sections["{ } Raw JSON"] = mo.json(dict(summary))
        return mo.vstack([header, mo.accordion(sections, multiple=True)], gap=0.8)

    # Plain AnnData summary.
    software = _software_badge(summary)
    header = _shape_row(summary, lead=[software] if software else None)
    sections = {
        "📋 Details": _anndata_panel(summary, show_shape=False),
        "{ } Raw JSON": mo.json(dict(summary)),
    }
    return mo.vstack([header, mo.accordion(sections, multiple=True)], gap=0.8)

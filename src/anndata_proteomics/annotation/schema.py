"""Pydantic models for the annotation TOML schema.

An annotation TOML attaches an external table to an AnnData/MuData axis. This first
iteration covers the ``obs`` (sample) axis only; the top-level shape leaves room for a
sibling ``[var]`` block (feature annotations) without changing the obs format.

Example::

    schema_version = "0.1"

    [obs]
    match_on  = "index"      # "index" => obs_names; else the name of an obs column to join on
    key_field = "raw_file"   # the field within each record that holds the join value

    [[obs.samples]]
    raw_file    = "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01"
    sample_name = "Condition_A_Sample_Alpha_01"
    condition   = "A"
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ObsAnnotation(_Strict):
    """Sample-annotation table joined onto ``obs``.

    Each ``samples`` record is free-form: every field other than ``key_field`` becomes an
    ``obs`` column. ``match_on`` selects what the record's ``key_field`` is matched against:
    ``"index"`` (the default) means ``obs_names``; any other value names an ``obs`` column.
    """

    match_on: str = "index"
    key_field: str = "raw_file"
    samples: list[dict[str, Any]] = Field(min_length=1)

    @model_validator(mode="after")
    def _records_carry_key_field(self) -> ObsAnnotation:
        for i, record in enumerate(self.samples):
            if self.key_field not in record:
                raise ValueError(
                    f"obs.samples[{i}] is missing key_field {self.key_field!r}; "
                    f"present fields: {sorted(record)}"
                )
        return self


class AnnotationSpec(_Strict):
    """Top-level annotation document."""

    schema_version: str
    obs: ObsAnnotation

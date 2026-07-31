"""Backend-independent FASTA workflow tests."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.annotation.validate_fasta import PeptideFastaMatchingConfig
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.workflows.fasta import (
    FeatureMappingAxes,
    PeptideLevelInput,
    PeptideLevelValidationWithReportedProteins,
    PeptideLevelWithReportedProteinsInput,
    PeptideValidationWorkflowResult,
    ProteinAnnotationInput,
    StoredProteinSearchParameters,
    build_mulink_feature_mapping,
    resolve_protein_annotation_input,
    validate_peptide_levels,
)

FASTA = ">sp|P1|ONE\nMPEPTIDEK\n>sp|P2|TWO\nMOTHERK\n"


def test_protein_digestion_is_resolved_from_explicit_physical_values() -> None:
    inputs = ProteinAnnotationInput(
        protein_groups=pd.Series(["P1"], index=["protein:P1"]),
        match_on="Protein_Group",
        search_parameters=StoredProteinSearchParameters(
            Parameters(
                enzyme="Lys-C",
                min_peptide_length=8,
                max_peptide_length=42,
            )
        ),
    )

    calculated = resolve_protein_annotation_input(inputs)

    assert calculated.config.fasta.cleavage.enzyme == "Lys-C"
    assert calculated.config.fasta.min_length == 8
    assert calculated.config.fasta.max_length == 42
    assert calculated.provenance.match_on == "Protein_Group"


def test_shared_validation_uses_only_typed_level_values() -> None:
    inputs = (
        PeptideLevelWithReportedProteinsInput(
            name="ion",
            sequences=pd.Series(["PEPTIDE"], index=["ion:peptide"]),
            reported_protein_field="Protein_Group",
            reported_proteins=pd.Series(["P1"], index=["ion:peptide"]),
        ),
        PeptideLevelInput(
            name="fragment",
            sequences=pd.Series(["OTHER"], index=["fragment:other"]),
        ),
    )

    calculated = validate_peptide_levels(
        inputs,
        FASTA,
        PeptideFastaMatchingConfig(),
    )

    ion = calculated.levels["ion"]
    assert isinstance(ion, PeptideLevelValidationWithReportedProteins)
    assert ion.reported_protein_field == "Protein_Group"
    assert ion.summary.loc["ion:peptide", "peptide_in_leading_protein"]
    assert pd.isna(
        calculated.levels["fragment"].summary.loc[
            "fragment:other",
            "leading_protein_in_fasta",
        ]
    )


def test_feature_mapping_is_calculated_from_explicit_axes() -> None:
    calculated = validate_peptide_levels(
        (
            PeptideLevelInput(
                name="ion",
                sequences=pd.Series(["PEPTIDE"], index=["ion:peptide"]),
            ),
        ),
        FASTA,
        PeptideFastaMatchingConfig(),
    )
    axes = FeatureMappingAxes(
        global_feature_names=pd.Index(["ion:peptide", "protein:P1"]),
        peptide_feature_names={"ion": pd.Index(["ion:peptide"])},
        protein_groups=pd.Series(["P1"], index=["protein:P1"]),
    )

    result = build_mulink_feature_mapping(axes, calculated)

    assert result.feature_mapping.mapping[0, 1] == 1
    assert result.target_rows.tolist() == [True, False]


def test_in_memory_read_and_write_functions_drive_the_workflow() -> None:
    source = {"ion": pd.Series(["PEPTIDE"], index=["ion:peptide"])}
    stored: dict[str, pd.DataFrame] = {}

    def read_inputs() -> tuple[PeptideLevelInput, ...]:
        return tuple(
            PeptideLevelInput(name=name, sequences=sequences.copy())
            for name, sequences in source.items()
        )

    def write_results(result: PeptideValidationWorkflowResult) -> None:
        stored.update(
            (name, validation.summary.copy()) for name, validation in result.levels.items()
        )

    calculated = validate_peptide_levels(
        read_inputs(),
        FASTA,
        PeptideFastaMatchingConfig(),
    )
    write_results(calculated)

    assert stored["ion"].loc["ion:peptide", "peptide_in_fasta"]
    assert source["ion"].tolist() == ["PEPTIDE"]

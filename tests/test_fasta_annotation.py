"""Tests for the FASTA → protein-annotation DataFrame builder."""

from __future__ import annotations

from io import StringIO

import pytest

from anndata_proteomics.fasta.annotation import (
    count_peptides,
    extract_gene_name,
    fasta_to_dataframe,
)
from anndata_proteomics.fasta.parser import iter_fasta

# Fixture string lifted from prolfquapp's get_annot_from_FASTA.R (.getSequences()).
# 7 forward records + 4 REV_ decoys = 11 records total.
PROLFQUAPP_FIXTURE = """\
>sp|A0A385XJL2|YGDT_ECOLI Protein YgdT OS=Escherichia coli (strain K12) OX=83333 GN=ygdT PE=4 SV=1
MLSTESWDNCEKPPLLFPFTALTCDETPVFSGSVLNLVAHSVDKYGIG
>sp|A5A615|YNCL_ECOLI Uncharacterized protein YncL OS=Escherichia coli (strain K12) OX=83333 GN=yncL PE=1 SV=1
MNVSSRTVVLINFFAAVGLFTLISMRFGWFI
>sp|P03018|UVRD_ECOLI DNA helicase II OS=Escherichia coli (strain K12) OX=83333 GN=uvrD PE=1 SV=1
MDVSYLLDSLNDKQREAVAAPRSNLLVLAGAGSGKTRVLVHRIAWLMSVENCSPYSIMAV
>sp|P04982|RBSD_ECOLI D-ribose pyranase OS=Escherichia coli (strain K12) OX=83333 GN=rbsD PE=1 SV=3
MKKGTVLNSDISSVISRLGHTDTLVVCDAGLPIPKSTTRIDMALTQGVPSFMQVLGVVTN
>sp|P04994|EX7L_ECOLI Exodeoxyribonuclease 7 large subunit OS=Escherichia coli (strain K12) OX=83333 GN=xseA PE=1 SV=2
MLPSQSPAIFTVSRLNQTVRLLLEHEMGQVWISGEISNFTQPASGHWYFTLKDDTAQVRC
>zz|Y-FGCZCont00001|  zz_FGCZCont0000_P61626_LYSC_HUMAN blastpHomologue_5.0e-107
MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRA
>zz|Y-FGCZCont00002|  zz_FGCZCont0001_P02534_K1M1_SHEEP blastpHomologue_0.0
SFNFCLPNLSFRSSCSSRPCVPSSCCGTTLPGACNIPANVGSCNWFCEGSFDGNEKETMQ
>REV_sp|Q13515|BFSP2_HUMAN Phakinin OS=Homo sapiens OX=9606 GN=BFSP2 PE=1 SV=1
GSEERDLLAHYSAVDKQLQCKRALLHAREQQQQEAEARIERLEAELRGVVAGLNQLEMDH
>REV_sp|Q14183|DOC2A_HUMAN Double C2-like domain-containing protein alpha OS=Homo sapiens OX=9606 GN=DOC2A PE=1 SV=5
ASSLAGAAPPLESTLTHWRELAADPQQLCDSWHKRAEGRAGPGLSVGGIFDNSKGIDYDW
>REV_tr|A0A075B6W8|A0A075B6W8_HUMAN T cell receptor alpha joining 17 (Fragment) OS=Homo sapiens OX=9606 GN=TRAJ17 PE=4 SV=1
PKVLVRTGGGFTLKNGAAKIX
>REV_sp|A0A385XJL2|YGDT_ECOLI Protein YgdT OS=Escherichia coli (strain K12) OX=83333 GN=ygdT PE=4 SV=1
GIGYKDVSHAVLNLVSGSFVPTEDCTLATFPFLLPPKECNDWSETSLM
"""


def test_parser_reads_wrapped_sequences():
    text = ">id1 head1\nACDE\nFGHI\n>id2 head2\nKLMN\n"
    records = list(iter_fasta(StringIO(text)))
    assert len(records) == 2
    assert records[0].header == "id1 head1"
    assert records[0].sequence == "ACDEFGHI"
    assert records[1].sequence == "KLMN"


def test_extract_gene_name_uniprot():
    header = "Protein YgdT OS=Escherichia coli (strain K12) OX=83333 GN=ygdT PE=4 SV=1"
    assert extract_gene_name(header) == "ygdT"


def test_extract_gene_name_non_uniprot_returns_empty():
    assert extract_gene_name("zz_FGCZCont0000_P61626_LYSC_HUMAN blastp") == ""


def test_count_peptides_matches_prolfquapp_docstring():
    # prolfquapp docstring example. Cleavage rule: K|R not followed by P
    # and not at end-of-string. Sequence "MKGLPRAKSHGSTGWGKRKRNKPK":
    #   cleavage ends (1-based): [2, 6, 8, 17, 18, 19, 20]
    #   segments → lengths: 2, 4, 2, 9, 1, 1, 1, 4
    # With min_length=5 only the length-9 peptide qualifies → 1.
    assert count_peptides("MKGLPRAKSHGSTGWGKRKRNKPK", min_length=5) == 1


def test_count_peptides_excludes_kr_before_proline():
    # KP and RP must not cleave. KR followed by anything else does cleave.
    # "AAAKPAAR" → K not cleaved (followed by P), R at the end (not cleaved).
    assert count_peptides("AAAKPAAR", min_length=1, max_length=100) == 1


def test_fasta_to_dataframe_columns_and_decoy_filter():
    df = fasta_to_dataframe(PROLFQUAPP_FIXTURE)
    # 4 REV_ decoys removed → 7 forward records
    assert len(df) == 7
    expected_cols = {
        "fasta.id",
        "fasta.header",
        "proteinname",
        "gene_name",
        "protein_length",
        "nr_peptides",
    }
    assert set(df.columns) == expected_cols

    first = df.iloc[0]
    assert first["fasta.id"] == "sp|A0A385XJL2|YGDT_ECOLI"
    assert "GN=ygdT" in first["fasta.header"]
    assert first["proteinname"] == "A0A385XJL2"
    assert first["gene_name"] == "ygdT"
    assert first["protein_length"] == 48


def test_fasta_to_dataframe_keeps_decoys_when_pattern_empty():
    df = fasta_to_dataframe(PROLFQUAPP_FIXTURE, decoy_pattern="")
    assert len(df) == 11


def test_fasta_to_dataframe_non_uniprot_proteinname_equals_id():
    df = fasta_to_dataframe(PROLFQUAPP_FIXTURE, is_uniprot=False)
    assert (df["proteinname"] == df["fasta.id"]).all()


def test_fasta_to_dataframe_gene_name_gated_on_match_count():
    # A 1-record UniProt-style fasta produces only one GN match → column omitted
    # to match prolfquapp's "added only if >1 matches" gating.
    text = ">sp|P0A6F5|GROEL_ECOLI Chaperone OS=Ecoli GN=groEL PE=1 SV=1\nAAAA\n"
    df = fasta_to_dataframe(text)
    assert "gene_name" not in df.columns


def test_fasta_to_dataframe_include_sequence():
    df = fasta_to_dataframe(PROLFQUAPP_FIXTURE, include_sequence=True)
    assert "sequence" in df.columns
    assert df.iloc[0]["sequence"].startswith("MLSTESW")


def test_fasta_to_dataframe_dedupes_on_fasta_id():
    # Same fasta.id repeated → first occurrence wins, second is dropped.
    text = (
        ">sp|P1|A first OS=x GN=a PE=1 SV=1\nAAAA\n"
        ">sp|P2|B second OS=x GN=b PE=1 SV=1\nCCCC\n"
        ">sp|P1|A third OS=x GN=c PE=1 SV=1\nDDDD\n"
    )
    df = fasta_to_dataframe(text, include_sequence=True)
    assert len(df) == 2
    assert df.loc[df["fasta.id"] == "sp|P1|A", "sequence"].iat[0] == "AAAA"


def test_fasta_to_dataframe_min_length_excludes_short_peptides():
    # Same record, two cutoffs → strictly fewer peptides for the higher floor.
    text = ">sp|P|X test GN=x PE=1 SV=1\nMKGLPRAKSHGSTGWGKRKRNKPK\n"
    df_lo = fasta_to_dataframe(text, min_length=2)
    df_hi = fasta_to_dataframe(text, min_length=5)
    assert df_lo["nr_peptides"].iat[0] > df_hi["nr_peptides"].iat[0]


@pytest.mark.parametrize(
    ("header", "expected_id", "expected_header"),
    [
        ("sp|P1|N foo", "sp|P1|N", "foo"),
        ("sp|P1|N;  foo bar", "sp|P1|N", "foo bar"),
        (">sp|P1|N foo", "sp|P1|N", "foo"),
    ],
)
def test_header_id_split_strips_leading_gt_and_trailing_semicolon(
    header, expected_id, expected_header
):
    text = f">{header}\nAAAA\n"
    df = fasta_to_dataframe(text)
    assert df["fasta.id"].iat[0] == expected_id
    assert df["fasta.header"].iat[0] == expected_header

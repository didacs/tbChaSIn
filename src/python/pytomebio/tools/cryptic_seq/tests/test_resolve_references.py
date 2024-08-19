from pathlib import Path
from typing import Dict

from pytomebio.tools.cryptic_seq.resolve_references import parse_reference
from pytomebio.tools.cryptic_seq.resolve_references import update_metasheet
from pytomebio.tools.cryptic_seq.utils import read_metasheet


def test_resolve_references(datadir: Path, references: Dict[str, str]) -> None:
    metasheet = datadir / "metasheet.xlsx"
    sample_df = read_metasheet(metasheet)
    parsed_references = {key: parse_reference(value) for key, value in references.items()}
    sample_df, unique_references = update_metasheet(
        sample_df, reference_map=parsed_references, root="/data/refs/"
    )
    refs = sample_df["reference"]
    assert all(refs[i] == "GRCh38.p14" for i in range(3))
    assert refs[3] == "hg19"
    assert unique_references == {
        "GRCh38.p14": "/data/refs/GRCh38/GRCh38.p14/",
        "hg19": "/data/refs/foo/bar/GRCh37/hg19/",
    }


def test_resolve_references_s3(datadir: Path, s3_references: Dict[str, str]) -> None:
    metasheet = datadir / "metasheet.xlsx"
    sample_df = read_metasheet(metasheet)
    parsed_references = {key: parse_reference(value) for key, value in s3_references.items()}
    sample_df, unique_references = update_metasheet(sample_df, reference_map=parsed_references)
    refs = sample_df["reference"]
    assert all(refs[i] == "GRCh38.p14" for i in range(3))
    assert refs[3] == "hg19"
    assert unique_references == {
        "GRCh38.p14": "s3://mock_bucket/GRCh38/GRCh38.p14/",
        "hg19": "s3://mock_bucket/foo/bar/GRCh37/hg19/",
    }

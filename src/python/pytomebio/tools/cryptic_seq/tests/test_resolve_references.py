from pathlib import Path
from typing import Dict

from pytomebio.tools.cryptic_seq.resolve_references import update_metasheet
from pytomebio.tools.cryptic_seq.utils import read_metasheet


def test_resolve_references(datadir: Path, references: Dict[str, str]) -> None:
    metasheet = datadir / "metasheet.xlsx"
    sample_df = read_metasheet(metasheet)
    sample_df, unique_references = update_metasheet(
        sample_df, reference_map=references, root="/data/refs/"
    )
    refs = sample_df["reference"]
    assert all(refs[i] == "GRCh38.p14" for i in range(3))
    assert refs[3] == "GRCh37"
    assert unique_references == {
        "/data/refs/GRCh38.p14/",
        "/data/refs/foo/bar/GRCh37/",
    }

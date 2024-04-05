from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from unittest.mock import MagicMock

import pandas as pd

from pytomebio.tools.cryptic_seq import create_metasheet_from_benchling as mod

CTB_ID = "123456789"
GENOMES = {"mouse": "mm10"}
DEFAULT_GENOME = "GRCh38/hg38"
ASSAY_QUERY_COLUMNS = [
    (name,)
    for name in (
        "group",
        "sample_name",
        "replicate",
        "integrase",
        "benchling_pl_id",
        "cell_type_or_clone",
        "species",
    )
]
ASSAY_QUERY_RESULT = [
    ("group1", "mickey", 1, "abc", "123", "brain", "mouse"),
    ("group2", "joe", 1, "abc", "456", "brain", "human"),
]
TRACKING_QUERY_COLUMNS = [
    (name,)
    for name in (
        "assay",
        "experiment_eln_id",
        "number_of_samples",
        "sequencing_run_id",
        "sequencing_project_name",
    )
]
TRACKING_QUERY_RESULT = ("foo", 123, 2, "abc123", "bar")


class MockCursor:
    def __init__(self, ctb_id: str) -> None:
        self._assay_query = mod.ASSAY_QUERY % {"ctb_id": ctb_id}
        self._tracking_query = mod.TRACKING_QUERY % {"ctb_id": ctb_id}
        self._last_query: str = None

    def execute(self, query: str, params: Dict[str, Any]) -> None:
        self._last_query = query % params

    @property
    def description(self) -> List[Tuple[Any, ...]]:
        if self._last_query == self._tracking_query:
            return TRACKING_QUERY_COLUMNS
        elif self._last_query == self._assay_query:
            return ASSAY_QUERY_COLUMNS
        else:
            raise ValueError(f"Unexpected query: {self._last_query}")

    def fetchone(self) -> Optional[Tuple[Any, ...]]:
        if self._last_query == self._tracking_query:
            return TRACKING_QUERY_RESULT
        else:
            raise ValueError(f"Unexpected query: {self._last_query}")

    def fetchall(self) -> List[Tuple[Any, ...]]:
        if self._last_query == self._assay_query:
            return ASSAY_QUERY_RESULT
        else:
            raise ValueError(f"Unexpected query: {self._last_query}")


def test_create_metasheet_from_benchling() -> None:
    mock_connection = MagicMock()
    mock_connection.cursor.return_value = MockCursor(ctb_id=CTB_ID)

    sample_df = mod.create_metasheet(
        mock_connection,
        ctb_id=CTB_ID,
        genomes=GENOMES,
        default_genome=DEFAULT_GENOME,
    )
    assert len(sample_df) == 2

    expected_df = pd.DataFrame(
        ASSAY_QUERY_RESULT, columns=[column[0] for column in ASSAY_QUERY_COLUMNS]
    )
    expected_df["genome_build"] = ("mm10", "GRCh38/hg38")
    assert sample_df.equals(expected_df)

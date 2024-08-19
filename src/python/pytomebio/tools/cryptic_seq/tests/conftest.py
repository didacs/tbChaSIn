from pathlib import Path
from typing import Dict

import pytest


@pytest.fixture
def datadir() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture
def attachment_sites() -> Dict[str, str]:
    ATT_SEQ = (
        "attP-GT:ACAGGCACAACCGTGGTTTGTCTGGTCAACCACCGCG:GT:CTCAGTGGTGTACGGTACAAACCCAGGATTAATAAGT"
    )
    return {"FRAG452": ATT_SEQ, "FRAG454": ATT_SEQ}


@pytest.fixture
def references() -> Dict[str, str]:
    return {
        "GRCh38/hg38": "GRCh38/GRCh38.p14/GRCh38.p14.fasta.gz",
        "GRCh37/hg37": "foo/bar/GRCh37/hg19/",
    }


@pytest.fixture
def s3_references() -> Dict[str, str]:
    return {
        "GRCh38/hg38": "s3://mock_bucket/GRCh38/GRCh38.p14/GRCh38.p14.fasta.gz",
        "GRCh37/hg37": "s3://mock_bucket/foo/bar/GRCh37/hg19/",
    }

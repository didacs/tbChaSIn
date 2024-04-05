from pathlib import Path
from typing import Dict

from pytomebio.tools.cryptic_seq.create_config_from_metasheet import create_groups
from pytomebio.tools.cryptic_seq.create_config_from_metasheet import update_metasheet
from pytomebio.tools.cryptic_seq.utils import read_metasheet


def test_create_config_from_metasheet(
    tmp_path: Path, datadir: Path, attachment_sites: Dict[str, str]
) -> None:
    sample_names = {
        "attB_NN_rep1_S3",
        "attB_NN_rep2_S5",
        "attB_NN_rep3_S9",
    }
    fastq_paths = {
        tmp_path / sample_name / f"reads_R{read}_001.fastq.gz"
        for sample_name in sample_names
        for read in (1, 2)
    }
    sample_names2 = {"attP_NN_rep3_S7"}
    fastq_paths2 = {
        tmp_path / sample_name / f"reads_R{read}_001.fastq.gz"
        for sample_name in sample_names2
        for read in (1, 2)
    }
    for path in fastq_paths | fastq_paths2:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    metasheet = datadir / "metasheet.txt"
    sample_df = read_metasheet(metasheet)
    samples_df = update_metasheet(
        sample_df,
        ref_dir=Path("references"),
        attachment_sites=attachment_sites,
        fastq_dir=tmp_path,
    )
    groups = create_groups(samples_df)
    assert len(groups) == 2
    groups_by_name = {group["name"]: group for group in groups}

    attB_NN = groups_by_name["attB_NN"]
    assert Path(attB_NN["ref_fasta"]) == Path("references/GRCh38.p14/GRCh38.p14.fasta")
    assert attB_NN["attachment_sites"] == [attachment_sites["FRAG452"]]
    assert len(attB_NN["samples"]) == 3
    assert set(sample["name"] for sample in attB_NN["samples"]) == sample_names
    assert set(sample["replicate"] for sample in attB_NN["samples"]) == {1, 2, 3}
    assert (
        set(Path(sample[col]) for sample in attB_NN["samples"] for col in ("fq1", "fq2"))
        == fastq_paths
    )

    attP_NN = groups_by_name["attP_NN"]
    assert Path(attP_NN["ref_fasta"]) == Path("references/GRCh37/GRCh37.fasta")
    assert attP_NN["attachment_sites"] == [attachment_sites["FRAG452"]]
    assert len(attP_NN["samples"]) == 1
    assert set(sample["name"] for sample in attP_NN["samples"]) == sample_names2
    assert set(sample["replicate"] for sample in attP_NN["samples"]) == {3}
    assert (
        set(Path(sample[col]) for sample in attP_NN["samples"] for col in ("fq1", "fq2"))
        == fastq_paths2
    )

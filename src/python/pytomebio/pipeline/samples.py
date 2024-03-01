from pathlib import Path
from typing import Any
from typing import Dict
from typing import List

import yaml
from attrs import define


@define
class Sample:
    """Stores sample data needed by the pytomebio tools.

    Attributes:
        name: The name of the sample
        group: The name of the group this sample is part of.
        replicate: the whole number (e.g. 1, 2, 3, ...) specifying the replicate for this sample
        fq1: The absolute path to the FASTQ for read 1 (R1)
        fq2: The absolute path to the FASTQ for read 2 (R2)
        ref_fasta: The absolute path to the reference FASTA, with accompanying BWA index files and
                   FASTA index
        attachment_sites: The list of attachment sites ("<name>:<left-seq>:<overhang>:<right-seq>")
    """

    name: str
    group: str
    replicate: int
    fq1: Path
    fq2: Path
    ref_fasta: Path
    attachment_sites: List[str]
    extra: Dict[str, Any]


def from_path(yml: Path) -> Dict[str, Sample]:
    """Returns dict from sample name to Sample for all the samples in a snakemake config.

    The snakemake config must be a yaml with a "settings" object in the top level.
    "settings" is a list of objects with group information:
        "name": the group name
        "ref_fasta": the absolute path to the reference used for every sample in this group
        "attachment_sites": the list of attachment sites used for every sample in this group
        "samples": a list of objects with specific information for each sample in the group
            "name": the name of each sample
            "replicate": the whole number (e.g. 1, 2, ...) specifying the replicate for this sample
            "fq1": The absolute path to the FASTQ for read 1 (R1)
            "fq2": The absolute path to the FASTQ for read 2 (R2)

    Args:
        yml: Path to yaml with snakemake config.
    """
    with yml.open("r") as reader:
        config: Dict[str, Any] = yaml.safe_load(reader)
    return _from_config(config=config)


def _from_config(config: Dict[str, Any]) -> Dict[str, Sample]:
    """Returns dict from sample name to Sample for all the samples in a snakemake config.

    The snakemake config must be a yaml with a "settings" object in the top level.
    "settings" is a list of objects with group information:
        "name": the group name
        "ref_fasta": the absolute path to the reference used for every sample in this group
        "attachment_sites": the list of attachment sites used for every sample in this group
        "samples": a list of objects with specific information for each sample in the group
            "name": the name of each sample
            "replicate": the whole number (e.g. 1, 2, ...) specifying the replicate for this sample
            "fq1": The absolute path to the FASTQ for read 1 (R1)
            "fq2": The absolute path to the FASTQ for read 2 (R2)

    Args:
        config: Dict object resulting from Path to yaml with snakemake config.
    """

    # Read in the samples from the config
    samples: List[Sample] = []
    for group in config["settings"]:
        for sample_args_dict in group["samples"]:
            # copy over values at the group level, and make sure to convert to the correct type
            sample_args: Dict[str, Any] = {
                "name": sample_args_dict["name"],
                "group": group["name"],
                "ref_fasta": Path(group["ref_fasta"]),
                "attachment_sites": group.get("attachment_sites", []),
                "fq1": Path(sample_args_dict["fq1"]),
                "fq2": Path(sample_args_dict["fq2"]),
                "replicate": int(sample_args_dict["replicate"]),
            }
            sample_args["extra"] = {
                key: value for key, value in sample_args_dict.items() if key not in sample_args
            }
            for key, value in group.items():
                if key == "samples":
                    continue
                sample_args["extra"][key] = value
            # create the sample and add it to the list
            sample = Sample(**sample_args)
            samples.append(sample)

    sample_dict: Dict[str, Sample] = {sample.name: sample for sample in samples}
    sample_names = [sample.name for sample in samples]
    assert len(sample_dict) == len(sample_names), "Sample names are not unique across ALL samples"

    return sample_dict

import json
from enum import Enum
from pathlib import Path
from typing import Any
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import yaml
from attrs import define

from pytomebio.core.sites import FindSitesMetric


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

    @classmethod
    def _from_config(cls, config: Dict[str, Any]) -> Dict[str, "Sample"]:
        """Returns dict from sample name to Sample for all the samples in a snakemake config.

        The snakemake config must be a yaml with a "settings" object in the top level.
        "settings" is a list of objects with group information:
            "name": the group name
            "ref_fasta": the absolute path to the reference used for every sample in this group
            "attachment_sites": the list of attachment sites used for every sample in this group
            "samples": a list of objects with specific information for each sample in the group
                "name": the name of each sample
                "replicate": the whole number (e.g. 1, 2, ...) specifying the replicate for this
                sample
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
                    if key != "samples":
                        sample_args["extra"][key] = value
                # create the sample and add it to the list
                sample = Sample(**sample_args)
                samples.append(sample)

        sample_dict: Dict[str, Sample] = {sample.name: sample for sample in samples}
        sample_names = [sample.name for sample in samples]
        assert len(sample_dict) == len(
            sample_names
        ), "Sample names are not unique across ALL samples"

        return sample_dict

    @classmethod
    def from_yml(cls, path: Path, ref_dir: Optional[Path] = None) -> Dict[str, "Sample"]:
        """Returns dict from sample name to Sample for all the samples in a snakemake config.

        The snakemake config must be a yaml with a "settings" object in the top level.
        "settings" is a list of objects with group information:
            "name": the group name
            "ref_fasta": the absolute path to the reference used for every sample in this group
            "attachment_sites": the list of attachment sites used for every sample in this group
            "samples": a list of objects with specific information for each sample in the group
                "name": the name of each sample
                "replicate": the whole number (e.g. 1, 2, ...) specifying the replicate for this
                sample
                "fq1": The absolute path to the FASTQ for read 1 (R1)
                "fq2": The absolute path to the FASTQ for read 2 (R2)

        Args:
            path: Path to yaml with snakemake config.
            ref_dir: Optional Path to the directory with reference FASTA files, if the reference
                paths in the YAML file are relative.
        """
        with path.open("r") as reader:
            config: Dict[str, Any] = yaml.safe_load(reader)
        return cls._from_config(config=config)


class FileGrouping(Enum):
    flat = "flat"
    group_name = "group_name"


@define
class SampleWithSites:
    """Stores information about a sample and its sites.

    Attributes:
        group: the sample group
        sample: the sample name
        replicate: the replicate number
        reference_filename: the name of the reference FASTA file
        metrics: the metrics associate with this sample
    """

    group: str
    name: str
    replicate: Optional[int] = None
    ref_fasta: Optional[Path] = None
    gene_set_gtf: Optional[Path] = None
    metrics_file: Optional[Path] = None
    metrics: Optional[List[FindSitesMetric]] = None

    REPLICATE_PATTERNS: ClassVar[List[str]] = ["_rep\\d+_", "-rep\\d-"]
    SAMPLE_NUMBER_PATTERNS: ClassVar[List[str]] = ["_S\\d+_", "_S\\d+$"]

    def get_metrics(self) -> List[FindSitesMetric]:
        if self.metrics is None and self.metrics_file:
            self.metrics = list(FindSitesMetric.read(self.metrics_file))
        if self.metrics is None:
            raise ValueError(f"No metrics found for {self.name}")
        return self.metrics

    @classmethod
    def _from_sample(
        cls, sample: Sample, grouping: FileGrouping, ref_dir: Optional[Path] = None
    ) -> "SampleWithSites":
        txt_name = f"{sample.name}.sites.txt"
        if grouping == FileGrouping.flat:
            txt = Path(txt_name)
        elif grouping == FileGrouping.group_name:
            txt = Path(sample.group) / sample.name / txt_name

        ref_fasta = sample.ref_fasta
        if ref_fasta and ref_dir and not ref_fasta.is_absolute():
            ref_fasta = ref_dir / ref_fasta

        gene_set_gtf = (
            Path(sample.extra["gene_set_gtf"]) if "gene_set_gtf" in sample.extra else None
        )
        if gene_set_gtf and ref_dir and not gene_set_gtf.is_absolute():
            gene_set_gtf = ref_dir / gene_set_gtf

        return SampleWithSites(
            group=sample.group,
            name=sample.name,
            replicate=sample.replicate,
            ref_fasta=ref_fasta,
            gene_set_gtf=gene_set_gtf,
            metrics_file=txt,
        )

    @classmethod
    def _from_dict(
        cls, sample_dict: Dict[str, Any], grouping: FileGrouping, ref_dir: Optional[Path] = None
    ) -> "SampleWithSites":
        txt_name = f"{sample_dict['name']}.sites.txt"
        if grouping == FileGrouping.flat:
            txt = Path(txt_name)
        elif grouping == FileGrouping.group_name:
            txt = Path(sample_dict["group"]) / sample_dict["name"] / txt_name

        ref_fasta = Path(sample_dict["ref_fasta"]) if "ref_fasta" in sample_dict else None
        if ref_fasta and ref_dir and not ref_fasta.is_absolute():
            ref_fasta = ref_dir / ref_fasta

        gene_set_gtf = Path(sample_dict["gene_set_gtf"]) if "gene_set_gtf" in sample_dict else None
        if gene_set_gtf and ref_dir and not gene_set_gtf.is_absolute():
            gene_set_gtf = ref_dir / gene_set_gtf

        return SampleWithSites(
            group=sample_dict["group"],
            name=sample_dict["name"],
            replicate=sample_dict["replicate"],
            ref_fasta=ref_fasta,
            gene_set_gtf=gene_set_gtf,
            metrics_file=txt,
        )

    @classmethod
    def from_yml(
        cls, path: Path, grouping: FileGrouping, ref_dir: Optional[Path] = None
    ) -> List["SampleWithSites"]:
        config_sample_dict: Dict[str, Sample] = Sample.from_yml(path=path)
        return [
            cls._from_sample(sample=s, grouping=grouping, ref_dir=ref_dir)
            for s in config_sample_dict.values()
        ]

    @classmethod
    def from_json(
        cls, path: Path, grouping: FileGrouping, ref_dir: Optional[Path] = None
    ) -> List["SampleWithSites"]:
        with open(path) as fd:
            return [cls._from_dict(d, grouping=grouping, ref_dir=ref_dir) for d in json.load(fd)]

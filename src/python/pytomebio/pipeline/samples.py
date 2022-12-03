import copy
from attrs import define
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List


@define
class Sample:
    name: str
    group: str
    fq_dir: Path
    ref_fasta: Path

    @property
    def fq1(self) -> Path:
        return self.fq_dir / f"{self.name}_R1_001.fastq.gz"

    @property
    def fq2(self) -> Path:
        return self.fq_dir / f"{self.name}_R2_001.fastq.gz"


def from_config(config: Dict[str, Any]) -> Dict[str, Sample]:
    """Returns the list of samples from a snakemake config.

    Todo: describe require structure.
    """

    # Read in the samples from the config
    samples: List[Sample] = []
    for group in config["settings"]:
        sample_names = group["samples"]
        args_dict = copy.deepcopy(group)
        del args_dict["samples"]
        args_dict["group"] = args_dict["name"]
        # convert keys that need to be type(Path)
        for key in ["fq_dir", "ref_fasta"]:
            args_dict[key] = Path(args_dict[key])
        for name in sample_names:
            # set sample name
            args_dict["name"] = name
            sample = Sample(**args_dict)
            samples.append(sample)

    sample_dict: Dict[str, Sample] = {sample.name: sample for sample in samples}
    sample_names = [sample.name for sample in samples]
    assert len(sample_dict) == len(sample_names), "Sample names are not unique across ALL samples"

    return sample_dict

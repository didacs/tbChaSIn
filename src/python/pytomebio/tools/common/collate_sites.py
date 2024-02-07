import statistics
from collections import Counter
from itertools import groupby
from pathlib import Path
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple

from attr import frozen
from fgpyo.util.metric import Metric

from pytomebio.pipeline.samples import Sample as ConfigSample
from pytomebio.pipeline.samples import from_path
from pytomebio.tools.change_seq.find_sites import FindSitesMetric


@frozen
class Sample:
    """Stores information about a sample inferred from the path

    Attributes:
        group: the sample group
        sample: the sample name
        replicate: the replicate number
        metrics: the metrics associate with this sample
    """

    group: str
    name: str
    replicate: Optional[int]
    metrics: List[FindSitesMetric]

    REPLICATE_PATTERNS: ClassVar[List[str]] = ["_rep\\d+_", "-rep\\d-"]
    SAMPLE_NUMBER_PATTERNS: ClassVar[List[str]] = ["_S\\d+_", "_S\\d+$"]

    @classmethod
    def build(cls, config_sample: ConfigSample) -> "Sample":
        txt = Path(f"{config_sample.group}/{config_sample.name}/{config_sample.name}.sites.txt")
        metrics: List[FindSitesMetric] = list(FindSitesMetric.read(txt))
        return Sample(
            group=config_sample.group,
            name=config_sample.name,
            replicate=config_sample.replicate,
            metrics=metrics,
        )


@frozen
class CollateSitesMetric(Metric["CollateSitesMetric"]):
    group: Optional[str]
    sample: str
    replicate: Optional[int]
    count_geq_1x: int
    count_geq_5x: int
    count_geq_10x: int
    count_geq_25x: int
    count_geq_50x: int
    count_geq_100x: int
    count_geq_1000x: int
    count_geq_5000x: int

    @classmethod
    def from_sample(cls, sample: Sample) -> "CollateSitesMetric":
        return CollateSitesMetric(
            group=sample.group,
            sample=sample.name,
            replicate=sample.replicate,
            count_geq_1x=sum(1 for m in sample.metrics if m.count >= 1),
            count_geq_5x=sum(1 for m in sample.metrics if m.count >= 5),
            count_geq_10x=sum(1 for m in sample.metrics if m.count >= 10),
            count_geq_25x=sum(1 for m in sample.metrics if m.count >= 25),
            count_geq_50x=sum(1 for m in sample.metrics if m.count >= 50),
            count_geq_100x=sum(1 for m in sample.metrics if m.count >= 100),
            count_geq_1000x=sum(1 for m in sample.metrics if m.count >= 1000),
            count_geq_5000x=sum(1 for m in sample.metrics if m.count >= 5000),
        )


def collate_sites(*, in_yml: Path, out_prefix: Path) -> None:
    """Collates detected integration sites across samples.

    The input config file must be the same config file used in the CHANGE-Seq and Cryptic-Seq
    Snakemake pipelines.  Please refer to the top-level README for the current format.

    The input sites file must contain a metrics file (from the `FindSitesMetric` class) and have
    the following path based on the values in the inpute config:
    `<sample-group>/<sample-name>/<sample-name>.sites.txt`

    The first table gives the # of sites with count >= X (column) for each sample (row).  This is
    output to `*.per_sample.txt`.

    The second table gives the # sites with count >= X (column) seen in exactly Y samples (row).
    This is output to `*.by_num_samples.all.txt`.

    The third table gives the same as table 2, but for only for a given set of replicates.  One
    output file is produced per replicate group.  These are output to:
    `*.replicates.<sample-name>.txt`.

    The fourth table unions all sites within each replicate group, taking the mean count for common
    sites, then computes same as table 2.  This is output to
    `*by_num_samples.unioned_replicates.txt`.

    The fifth table outputs useful statistics across all detected sites.  This is output to:
    `*collated.txt`.

    Args:
        in_yml: the input configuration YAML used in the CHANGE-Seq or Cryptic-Seq Snakemake
                pipeline.
        out_prefix: the output prefix for all output files.
    """
    # Load in the samples using the pipeline config
    config_sample_dict: Dict[str, ConfigSample] = from_path(yml=in_yml)
    samples: List[Sample] = [Sample.build(config_sample=s) for s in config_sample_dict.values()]

    # Table 1
    CollateSitesMetric.write(
        Path(f"{out_prefix}.per_sample.txt"),
        *[CollateSitesMetric.from_sample(sample=sample) for sample in samples],
    )

    # Table 2
    CollateSitesMetric.write(
        Path(f"{out_prefix}.by_num_samples.all.txt"), *compute_metrics(samples=samples)
    )

    # Table 3
    for (group, name), replicates in groupby(samples, lambda sample: (sample.group, sample.name)):
        CollateSitesMetric.write(
            Path(f"{out_prefix}.replicates.{group}.{name}.txt"),
            *compute_metrics(samples=list(replicates)),
        )

    # Table 4
    unioned_replicates: List[Sample] = [
        union_replicates(replicates=list(replicates))
        for _, replicates in groupby(samples, lambda sample: (sample.group, sample.name))
    ]
    CollateSitesMetric.write(
        Path(f"{out_prefix}.by_num_samples.unioned_replicates.txt"),
        *compute_metrics(samples=unioned_replicates),
    )

    # Table 5
    site_counter: Dict[Tuple[str, int], List[int]] = {}
    site_to_samples: Dict[Tuple[str, int], List[str]] = {}
    for sample in samples:
        for metric in sample.metrics:
            key = (metric.reference_name, metric.position)
            if key not in site_counter:
                site_counter[key] = [metric.count]
                site_to_samples[key] = [sample.name]
            else:
                site_counter[key].append(metric.count)
                site_to_samples[key].append(sample.name)
    means: Dict[Tuple[str, int], float] = {
        ref_and_pos: (sum(counts) / float(len(counts)))
        for ref_and_pos, counts in site_counter.items()
    }
    with Path(f"{out_prefix}.collated.txt").open("w") as writer:
        writer.write(
            "reference_name\tposition\tnum_samples\tmean_count\t"
            "stddev_count\tmin_count\tmax_count\tsamples\n"
        )
        for ref_and_pos, mean_count in sorted(means.items(), key=lambda tup: -tup[1]):
            counts: List[int] = site_counter[ref_and_pos]
            reference_name, position = ref_and_pos
            sample_count = len(counts)
            stddev_count = statistics.stdev(counts) if sample_count > 1 else 0
            min_count = min(counts)
            max_count = max(counts)
            writer.write(
                f"{reference_name}\t{position}\t{sample_count}\t{mean_count:.1f}"
                f"\t{stddev_count:.1f}\t{min_count}\t{max_count}\t"
            )
            writer.write(",".join(site_to_samples[ref_and_pos]))
            writer.write("\n")


def union_replicates(replicates: List[Sample]) -> Sample:
    sample_counter: Counter = Counter()
    total_counter: Counter = Counter()
    for replicate in replicates:
        for metric in replicate.metrics:
            sample_counter[(metric.reference_name, metric.position)] += 1
            total_counter[(metric.reference_name, metric.position)] += metric.count

    metrics: List[FindSitesMetric] = []
    for loci in sample_counter:
        reference_name, position = loci
        sample_count = sample_counter[loci]
        total_count = total_counter[loci]
        mean_count = total_count // sample_count
        metric = FindSitesMetric(
            reference_name=reference_name,
            position=position,
            count=mean_count,
            attachment_site="NA",
            positive_strand=True,
            left_to_count=Counter(),
            right_to_count=Counter(),
            both_to_count=Counter(),
        )
        metrics.append(metric)
    sample: Sample = Sample(
        group=replicates[0].group, name=replicates[0].name, replicate=None, metrics=metrics
    )
    return sample


def compute_metrics(samples: List[Sample]) -> List[CollateSitesMetric]:
    min_count_to_counter: Dict[int, Counter] = {
        min_count: compute_sample_counter(samples=samples, min_count=min_count)
        for min_count in [1, 5, 10, 25, 50, 100, 1000, 5000]
    }

    metrics: List[CollateSitesMetric] = []
    for num_samples in range(1, len(samples) + 1):
        metric = CollateSitesMetric(
            group=None,
            sample=f"{num_samples}",
            replicate=None,
            count_geq_1x=min_count_to_counter[1][num_samples],
            count_geq_5x=min_count_to_counter[5][num_samples],
            count_geq_10x=min_count_to_counter[10][num_samples],
            count_geq_25x=min_count_to_counter[25][num_samples],
            count_geq_50x=min_count_to_counter[50][num_samples],
            count_geq_100x=min_count_to_counter[100][num_samples],
            count_geq_1000x=min_count_to_counter[1000][num_samples],
            count_geq_5000x=min_count_to_counter[5000][num_samples],
        )
        metrics.append(metric)
    return metrics


def compute_sample_counter(samples: List[Sample], min_count: int) -> Counter:
    # get the number of samples seen per site with at least min_count coverage
    site_counter: Counter = Counter()
    for sample in samples:
        for metric in sample.metrics:
            if metric.count >= min_count:
                site_counter[(metric.reference_name, metric.position)] += 1

    # return the number of sites seen in exactly X samples with at least min_count coverage
    sample_counter: Counter = Counter()
    for num_samples in range(1, len(samples) + 1):
        sample_counter[num_samples] = 0
    for _, num_samples in site_counter.items():
        sample_counter[num_samples] += 1
    return sample_counter

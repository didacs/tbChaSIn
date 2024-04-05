import glob
import json
import sys
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

import pandas as pd
import yaml

from pytomebio.tools.cryptic_seq.utils import read_metasheet
from pytomebio.tools.cryptic_seq.utils import resolve_local_path

DEFAULT_FASTQ_NAME_PATTERN = "**/{sample_name}*/*_R{read}_*.gz"


def _add_ref_fasta(row: pd.Series, *, ref_dir: Path) -> pd.Series:
    reference = row["reference"]
    row["ref_fasta"] = str(ref_dir / reference / f"{reference}.fasta")
    return row


def _add_attachment_sites(
    row: pd.Series,
    *,
    attachment_sites: Dict[str, str],
    default_attachment_site: Optional[str] = None,
) -> pd.Series:
    """
    If the `attachment_site` column is not specified for `row`, then it is set it to the
    concatenation of the attachment sites associated with the IDs in the `benchling_pl_id` column.

    Args:
        row: a row from the metasheet.
        attachment_sites: a mapping of attachment site IDs to sequences.
    """
    attachment_site = row.get("attachment_site")
    row_attachment_sites = None
    if pd.isna(attachment_site):
        ids = row["benchling_pl_id"].split(";")
        row_attachment_sites = [attachment_sites[id] for id in ids if id in attachment_sites]
    if row_attachment_sites:
        row["attachment_site"] = ";".join(row_attachment_sites)
    elif default_attachment_site is not None:
        row["attachment_site"] = default_attachment_site
    else:
        raise ValueError(f"No valid attachment site(s) specified for sample {row.sample_name}")

    return row


def _resolve_glob(pattern: str, *, root: Optional[Path] = None) -> str:
    """
    Resolves a path that may contain a glob expression. If `root` is specified, then `pattern` is
    expected to be relative to the root directory, otherwise it is expected to resolve to an
    absolute path.
    """
    if root is not None:
        pattern = str(Path(root) / pattern)

    paths = glob.glob(pattern, recursive=True)
    if len(paths) != 1:
        raise ValueError(f"Expected exactly one file matching {pattern}")

    return str(resolve_local_path(Path(paths[0])))


def _add_fastq_paths(
    row: pd.Series,
    *,
    fastq_dir: Optional[Path] = None,
    fastq_name_pattern: str = DEFAULT_FASTQ_NAME_PATTERN,
    paired_end: bool = True,
) -> pd.Series:
    """Adds 'fq1' and 'fq2' columns to the specified `row` with paths to the FASTQ files specified
    by the combination of the given `fastq_dir` and `fastq_name_pattern`. The format string can
    include `{column}` to substitute the value from the column in `row`. If the file name contains
    the read number (1 or 2), then the pattern must contain `{read}`.

    Args:
        row: a row from the metasheet.
        fastq_dir: the path to the base directory where FASTQ files are stored. If not specified,
                   then `fastq_name_format` is expected to be an absolute path.
        fastq_name_pattern: the format of the FASTQ file names, relative to `fastq_dir`.
        paired_end: whether the FASTQ files are paired-end or single-end.
    """
    fq1_name = row.get("fq1")
    if fq1_name is None:
        fq1_name = fastq_name_pattern.format(read=1, **row)

    fq1_path = _resolve_glob(fq1_name, root=fastq_dir)

    if not paired_end:
        fq2_path = None
    else:
        fq2_name = row.get("fq2")
        if fq2_name is None:
            fq2_name = fastq_name_pattern.format(read=2, **row)

        fq2_path = _resolve_glob(fq2_name, root=fastq_dir)
        if fq2_path == fq1_path:
            raise ValueError(f"FASTQ files {fq1_path} and {fq2_path} are the same")

    row["fq1"] = fq1_path
    row["fq2"] = fq2_path
    return row


def update_metasheet(
    sample_df: pd.DataFrame,
    *,
    ref_dir: Path,
    attachment_sites: Dict[str, str],
    default_attachment_site: Optional[str] = None,
    fastq_dir: Optional[Path] = None,
    fastq_name_pattern: str = DEFAULT_FASTQ_NAME_PATTERN,
) -> pd.DataFrame:
    """
    Parses a metasheet excel file and returns a pandas dataframe. Creates/fills the
    `attachment_site`, `fq1`, and `fq2` columns if they are missing/empty.

    Args:
        sample_df: the sample DataFrames.
        ref_dir: the path to the directory containing the references.
        attachment_sites: mapping of attachment site names to sequences.
        default_attachment_site: the default value for the `attachment_site` column.
        fastq_dir: the path to the base directory where FASTQ files are stored. If not specified,
            then `fastq_name_pattern` is expected to be an absolute path.
        fastq_name_pattern: the format of the FASTQ file names, relative to `fastq_dir. The format
            string can contain '{column}' to substitute the value from the column in the metasheet.
             If the file name contains the read number (1 or 2), then the pattern must contain
             '{read}'.
    """
    # add the ref_dir column
    sample_df = sample_df.apply(
        _add_ref_fasta,
        ref_dir=ref_dir,
        axis=1,
    )
    # add the attachment sites based on the benchling IDs
    sample_df = sample_df.apply(
        _add_attachment_sites,
        attachment_sites=attachment_sites,
        default_attachment_site=default_attachment_site,
        axis=1,
    )
    # add the fastq paths based on the specified base directory and name pattern
    sample_df = sample_df.apply(
        _add_fastq_paths, fastq_dir=fastq_dir, fastq_name_pattern=fastq_name_pattern, axis=1
    )
    # make sure there are no missing values
    if sample_df.isna().any().any():
        raise ValueError("Metasheet has missing values")
    return sample_df


def create_groups(sample_df: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Creates the settings groups (the value of the `settings` field in the config file) from a
    pandas dataframe of sample information. Returns a list of groups, where a group is determined
    by the value of the `group` column and the value is a dict of group-level settings.

    Args:
        sample_df: a pandas dataframe of sample information.
    """

    def set_or_compare(group: Dict[str, Any], key: str, value: Any) -> None:
        if key not in group:
            group[key] = value
        elif group[key] != value:
            raise ValueError(
                f"Samples with same group had different values for {key}: {group[key]} != {value}"
            )

    groups: Dict[str, dict] = {}
    for _, row in sample_df.iterrows():
        group_name = row["group"]
        group = groups.setdefault(group_name, {})
        # Set/check the group-level metadata
        set_or_compare(group, "name", group_name)
        # TODO: ref_fasta is going to cause trouble when we try to run this on AWS
        # We'll need to have a mapping of reference name to file instead
        set_or_compare(group, "ref_fasta", row["ref_fasta"])
        set_or_compare(group, "attachment_sites", set(row["attachment_site"].split(";")))
        # Add the sample-level metadata
        samples = group.setdefault("samples", {})
        # Check for duplicate sample names
        if row.sample_name in samples:
            raise ValueError(f"Duplicate sample name {row.sample_name}")
        samples[row.sample_name] = {
            "name": row["sample_name"],
            "replicate": row["replicate"],
            "fq1": row["fq1"],
            "fq2": row["fq2"],
        }

    def make_serializable(group: Dict[str, Any]) -> Dict[str, Any]:
        group["attachment_sites"] = list(group["attachment_sites"])
        group["samples"] = list(group["samples"].values())
        return group

    return [make_serializable(group) for group in groups.values()]


def create_config_from_metasheet(
    *,
    metasheet: Path,
    ref_dir: Path = Path("."),
    attachment_sites: Optional[Path] = None,
    default_attachment_site: Optional[str] = None,
    fastq_dir: Optional[Path] = None,
    fastq_name_pattern: str = DEFAULT_FASTQ_NAME_PATTERN,
    output_file: Optional[Path] = None,
    groups_file: Optional[Path] = None,
) -> None:
    """Creates the Cryptic-seq sample config file from a metasheet.

    Args:
        metasheet: the path to the metasheet.
        ref_dir: the path containing the references. Each reference must be a folder containing the
            FASTA file and BWA index files.
        attachment_sites: the path to the JSON file with a mapping of attachment site names to
            sequences.
        default_attachment_site: the default value for the `attachment_site` column.
        fastq_dir: the path to the base directory where FASTQ files are stored. If not specified,
            then --fastq-name-pattern is expected to be an absolute path.
        fastq_name_pattern: the format of the FASTQ file names, relative to --fastq-dir. The format
            string can contain '{column}' to substitute the value from the column in the
            metasheet. If the file name contains the read number (1 or 2), then the pattern must
            contain '{read}'.
        output_file: the path to the output config file. If not specified, then the config file is
            written to stdout.
        groups_file: the path to a file where the group names used in the config file will be
            written, one per line.
    """
    sample_df = read_metasheet(metasheet)

    attachment_sites_dict = {}
    if attachment_sites is not None:
        with open(attachment_sites) as fd:
            attachment_sites_dict = json.load(fd)

    sample_df = update_metasheet(
        sample_df,
        ref_dir=ref_dir,
        attachment_sites=attachment_sites_dict,
        default_attachment_site=default_attachment_site,
        fastq_dir=fastq_dir,
        fastq_name_pattern=fastq_name_pattern,
    )

    groups = create_groups(sample_df)

    # make each group serializable: convert attachment sites and samples to lists
    config = {"settings": groups}

    # write the config to a file, or to stdout if `output_file` is `None`
    if output_file is not None:
        with open(output_file, "w") as fd:
            yaml.dump(config, fd, indent=2, sort_keys=False)
    else:
        yaml.dump(config, sys.stdout, indent=2, sort_keys=False)

    if groups_file is not None:
        with open(groups_file, "w") as fd:
            fd.write("\n".join(group["name"] for group in groups))
            fd.write("\n")

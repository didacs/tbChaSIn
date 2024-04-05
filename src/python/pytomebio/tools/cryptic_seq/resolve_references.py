import json
import sys
from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Set
from typing import Tuple
from urllib.parse import urlparse

import pandas as pd

from pytomebio.tools.cryptic_seq.utils import read_metasheet
from pytomebio.tools.cryptic_seq.utils import resolve_local_path


def _check_url(url: str, *, require_dir: bool = False) -> Tuple[str, bool]:
    """
    Parses the given string as a URL. If it is a file (no scheme, or a scheme of 'file'),
    then returns `(path, False)` where `path` is the path component of the URL, otherwise
    returns `(url, True)`.
    """
    parsed = urlparse(url)
    if parsed.scheme in ["", "file"]:
        if parsed.netloc != "":
            raise ValueError(
                f"Invalid file URL {url}; file URLs must be of the form 'file:///path/to/file'"
            )
        return parsed.path, False
    elif not require_dir or url.endswith("/"):
        return url, True
    else:
        raise ValueError(f"URL {url} must end with a slash")


def _resolve_reference(path: str, *, root: Optional[str] = None) -> Tuple[str, str]:
    """
    Given a path and and optional root folder/URL, resolves the path to an absolute path. If `root`
    is `None` then `path` may be relative or absolute. Returns `(abs_path, name)` where `abs_path`
    is the absolute path and `name` is the file name. Raises an exception if an absolute path
    cannot be resolved.
    """
    if root is None:
        abs_path = str(resolve_local_path(Path(path)))
    else:
        root, root_is_url = _check_url(root, require_dir=True)
        if root_is_url:
            abs_path = f"{root}{path}"
        else:
            abs_path = str(resolve_local_path(Path(root) / path))

    last_slash_idx = abs_path.rfind("/")
    assert last_slash_idx >= 0
    return f"{abs_path}/", abs_path[last_slash_idx + 1 :]


def _set_reference(
    row: pd.Series,
    *,
    unique: Dict[str, str],
    reference_map: Dict[str, str],
    root: Optional[str] = None,
    default: Optional[str] = None,
) -> pd.Series:
    """
    If the `reference` column is present, it may be a key in the `reference_map`, a folder name,
    a relative path, or an absolute path. We attempt to resolve it to an absolute path. The
    column is then updated to just the folder name, and the absolute path is added to `unique`. If
    the `reference` column is not present, the `genome_build` column is looked up in
    `reference_map`. If it is not found, then the default reference is used. A `ValueError` is
    raised if the reference cannot be resolved.
    """
    reference = row.get("reference")
    genome = row.get("genome_build")
    print(f"reference: {reference} genome: {genome}")
    if not pd.isna(reference) and reference in reference_map:
        reference = reference_map[reference]
    elif not pd.isna(genome) and genome in reference_map:
        reference = reference_map[genome]
    elif default is not None:
        reference = default
    else:
        raise ValueError(f"No valid reference specified for sample {row.sample_name}")

    abs_path, name = _resolve_reference(reference, root=root)
    existing = unique.get(name)
    if existing is None:
        unique[name] = abs_path
    elif existing != abs_path:
        raise ValueError(f"Multiple references with name {name}: {unique.get(name)}, {abs_path}")

    row["reference"] = name
    return row


def update_metasheet(
    sample_df: pd.DataFrame,
    root: Optional[str] = None,
    reference_map: Dict[str, str] = None,
    default: Optional[str] = None,
) -> Tuple[pd.DataFrame, Set[str]]:
    unique_references: Dict[str, str] = {}
    sample_df = sample_df.apply(
        _set_reference,
        unique=unique_references,
        reference_map=reference_map,
        root=root,
        default=default,
        axis=1,
    )
    return sample_df, set(unique_references.values())


def resolve_references(
    *,
    metasheet: Path,
    root: Optional[str] = None,
    mapping: Optional[Path] = None,
    default: Optional[str] = None,
    output_metasheet: Optional[Path] = None,
    output_paths: Optional[Path] = None,
) -> None:
    """Adds/updates the 'genome_build' column in the specified metasheet, and outputs a new
    metasheet and a list of unique reference paths.

    Args:
        metasheet: the path to the input metasheet.
        root: the path to the root folder where references are stored. Each reference is a folder
            that contains the FASTA file and BWA index files.
        mapping: the path to the JSON file with a mapping of genomes to files.
        default: the path to the default reference fasta file.
        output_metasheet: the path to write the updated metasheet. Defaults to stdout.
        output_paths: the path to write a file with the unique reference paths (one per line).
            Defaults to './paths.txt'.
    """

    reference_map = {}
    if mapping is not None:
        with open(mapping) as fd:
            reference_map = json.load(fd)

    sample_df = read_metasheet(metasheet)
    sample_df, unique_references = update_metasheet(
        sample_df,
        reference_map=reference_map,
        root=root,
        default=default,
    )

    if output_metasheet is not None:
        sample_df.to_csv(output_metasheet or sys.stdout, sep="\t", index=False)

    with open(output_paths or "paths.txt", "w") as fd:
        fd.write("\n".join(sorted(unique_references)))
        fd.write("\n")

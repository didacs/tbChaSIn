import json
import sys
from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Tuple
from typing import Union
from urllib.parse import urlparse

import pandas as pd

from pytomebio.tools.cryptic_seq.utils import read_metasheet
from pytomebio.tools.cryptic_seq.utils import resolve_local_path

FASTA_EXT = ".fasta.gz"


def _as_folder(path: Union[str, Path]) -> str:
    path_str = str(path)
    return path_str if path_str.endswith("/") else f"{path_str}/"


def _check_folder_url(url: str, *, require_dir: bool = False) -> Tuple[str, bool]:
    """
    Parses the given string as a folder URL. If it is a local path (no scheme, or a scheme of
    'file'), then returns `(path, False)` where `path` is the path component of the URL, otherwise
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


def resolve_reference(folder_path: str, root: Optional[str] = None) -> str:
    """
    Given a root folder/URL, assumes `folder_path` is relative and resolves it to an absolute
    path. If `root` is `None` then `folder_path` may be relative or absolute and is resolved
    relative to the working directory. Returns the absolute path to the reference folder. Raises an
    exception if an absolute path cannot be resolved.
    """
    if root is None:
        return _as_folder(resolve_local_path(Path(folder_path)))
    else:
        root, root_is_url = _check_folder_url(root, require_dir=True)
        if root_is_url:
            return f"{root}{folder_path}"
        else:
            return _as_folder(resolve_local_path(Path(root) / folder_path))


def get_reference_id(folder_path: str) -> str:
    """
    Returns the ID of this reference, which is the last component of the folder path.
    """
    assert folder_path.endswith("/")
    last_slash_idx = folder_path.rfind("/", 0, -1)
    start = (last_slash_idx + 1) if last_slash_idx >= 0 else 0
    return folder_path[start:-1]


def parse_reference(path: str) -> str:
    """
    Parses the given string as a path/URI to either a folder or a .fasta.gz file. If it is a
    file, then checks that the file name without the .fasta.gz extension is the same as the
    last component of the folder path.
    """
    if path.endswith("/"):
        return path
    elif path.endswith(FASTA_EXT):
        # it is a file, split into folder and prefix
        last_slash_idx = path.rfind("/")
        if last_slash_idx >= 0:
            folder = path[: last_slash_idx + 1]
            prefix = path[last_slash_idx + 1 : -len(FASTA_EXT)]
            if not folder.endswith(f"{prefix}/"):
                raise ValueError(
                    f"Reference folder {folder} does not match the file prefix {prefix}"
                )
            return folder
        else:
            prefix = path[: -len(FASTA_EXT)]
            return f"{prefix}/"
    else:
        raise ValueError(f"Reference {path} is not a folder or a file with extension {FASTA_EXT}")


def _set_or_check_reference(
    reference_path: str,
    resolved: Dict[str, str],
    root: Optional[str] = None,
) -> str:
    resolved_reference_path = resolve_reference(reference_path, root)
    reference_id = get_reference_id(resolved_reference_path)
    existing_reference_path = resolved.get(reference_id)
    if existing_reference_path is None:
        resolved[reference_id] = resolved_reference_path
    elif existing_reference_path != resolved_reference_path:
        raise ValueError(
            "Multiple different references with name "
            f"{reference_id}: {existing_reference_path}, {resolved_reference_path}"
        )
    return reference_id


def _set_reference(
    row: pd.Series,
    *,
    resolved: Dict[str, str],
    reference_map: Dict[str, str],
    root: Optional[str] = None,
    default: Optional[str] = None,
) -> pd.Series:
    """
    If the `genome_build` column is present, the value is looked up in `reference_map`. If it is
    not found, then the default reference is used. The default may be a name in `reference_map` or
    a path. A `ValueError` is raised if the reference cannot be resolved.

    The reference ID is the last component of the reference folder path. The `reference` column is
    added/updated with the reference ID, and the mapping from reference ID to reference is added to
    the `resolved` dict if it is not already present.
    """
    genome = row.get("genome_build")
    has_genome = not pd.isna(genome)
    if has_genome and genome in reference_map:
        reference_path = reference_map[genome]
    elif default is not None:
        if default in reference_map:
            reference_path = reference_map[default]
        else:
            reference_path = default
    elif has_genome:
        reference_path = genome
    else:
        raise ValueError(
            "Either 'genome_build' must be specified in the metasheet, or a default reference "
            "must be provided"
        )

    row["reference"] = _set_or_check_reference(reference_path, resolved, root=root)
    return row


def update_metasheet(
    sample_df: pd.DataFrame,
    reference_map: Dict[str, str],
    root: Optional[str] = None,
    default: Optional[str] = None,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    resolved_references: Dict[str, str] = {}
    sample_df = sample_df.apply(
        _set_reference,
        resolved=resolved_references,
        reference_map=reference_map,
        root=root,
        default=default,
        axis=1,
    )
    return sample_df, resolved_references


def set_annotation_reference(
    resolved: Dict[str, str],
    reference_map: Dict[str, str],
    root: Optional[str] = None,
    annotation: Optional[str] = None,
    default: Optional[str] = None,
) -> str:
    ref = annotation or default
    if ref in reference_map:
        ref = reference_map[ref]

    if ref:
        return _set_or_check_reference(ref, resolved, root=root)
    elif len(resolved) == 1:
        return list(resolved.keys())[0]
    else:
        raise ValueError("Could not determine the annotation reference.")


def resolve_references(
    *,
    metasheet: Path,
    root: Optional[str] = None,
    mapping: Optional[Path] = None,
    default: Optional[str] = None,
    annotation: Optional[str] = None,
    output_metasheet: Optional[Path] = None,
    output_references: Optional[Path] = None,
    output_annotation: Optional[str] = None,
) -> None:
    """Adds/updates the 'reference' column in the specified metasheet, and outputs a new
    metasheet and a JSON file mapping each genome to the full path of its reference.

    Args:
        metasheet: the path to the input metasheet.
        root: the path to the root folder where references are stored. Each reference is a folder
            that contains the FASTA file and BWA index files.
        mapping: the path to the JSON file with a mapping of genomes to files.
        default: the name of, or path to, the default reference.
        annotation: the name of, or path to, the annotation reference.
        output_metasheet: the path to write the updated metasheet. Defaults to stdout.
        output_references: the path to write a JSON file with the mapping of genome to reference
            paths. Defaults to './resolved_references.json'.
        output_annotation: the path to write a text file with the name of the annotation reference.
            Defaults to './annotation_reference.txt'.
    """

    reference_map = {}
    if mapping is not None:
        with open(mapping) as fd:
            reference_map = {key: parse_reference(ref) for key, ref in json.load(fd).items()}

    sample_df = read_metasheet(metasheet)
    sample_df, resolved_references = update_metasheet(
        sample_df,
        reference_map=reference_map,
        root=root,
        default=default,
    )

    annotation_reference_id = set_annotation_reference(
        resolved_references,
        reference_map=reference_map,
        root=root,
        annotation=annotation,
        default=default,
    )

    if output_metasheet is not None:
        sample_df.to_csv(output_metasheet or sys.stdout, sep="\t", index=False)

    with open(output_references or "resoved_references.json", "w") as fd:
        json.dump(resolved_references, fd, indent=2)

    with open(output_annotation or "annotation_reference.txt", "w") as fd:
        fd.write(annotation_reference_id)

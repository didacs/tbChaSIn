from pathlib import Path

import pandas as pd


def read_metasheet(path: Path) -> pd.DataFrame:
    # read in the metasheet
    if path.suffix.startswith(".xls"):
        return pd.read_excel(path)
    else:
        return pd.read_csv(path, sep="\t", index_col=False)


def resolve_local_path(path: Path) -> Path:
    path = path.resolve()
    if not path.is_absolute():
        raise ValueError(f"{path} is not an absolute path")
    return path

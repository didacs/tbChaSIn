from pathlib import Path


def touch_path(path: Path) -> Path:
    """Touch the given path, creating parent directories if necessary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w"):
        pass
    return path

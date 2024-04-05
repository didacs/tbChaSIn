import json
import logging
import os
import sys
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Optional

import pandas as pd
import psycopg as pg
from psycopg.connection import Connection

LOG = logging.getLogger(__name__)

ENV_WAREHOUSE_USERNAME = "WAREHOUSE_USERNAME"
ENV_WAREHOUSE_PASSWORD = "WAREHOUSE_PASSWORD"
ENV_WAREHOUSE_HOST = "WAREHOUSE_HOST"
ENV_WAREHOUSE_PORT = "WAREHOUSE_PORT"
DEFAULT_WAREHOUSE_PORT = 5432
ENV_WAREHOUSE_DBNAME = "WAREHOUSE_DBNAME"
DEFAULT_WAREHOUSE_DBNAME = "warehouse"
ENV_WAREHOUSE_SSLMODE = "WAREHOUSE_SSLMODE"
DEFAULT_WAREHOUSE_SSLMODE = "verify-ca"

# Queries the genomic_assays_metadata table on Benchling. Should return one row for each sample
# with the specified ctb_id.
ASSAY_QUERY = """
select
    ga.group,
    ga.sample_name,
    ga.replicate,
    ga.protein_integrase as integrase,
    ga.donor as benchling_pl_id,
    ga.cell_type_or_clone,
    ga.species
from genomic_assays_metadata$raw as ga
where ga.archived$ = false and ga.ctb_id = '%(ctb_id)s'
"""

# Queries the custom_tracking table on Benchling. Should return one row for the specified ctb_id.
TRACKING_QUERY = """
select
    ct.assay,
    ct.experiment_eln_id,
    ct.number_of_samples,
    ct.sequencing_run_id,
    ct.sequencing_project_name
from custom_tracking$raw as ct
where ct.archived$ = false and ct.file_registry_id$ = '%(ctb_id)s'
"""


def warehouse_connect(
    username: Optional[str] = None,
    password: Optional[str] = None,
    host: Optional[str] = None,
    port: Optional[int] = None,
    dbname: Optional[str] = None,
    sslmode: Optional[str] = None,
) -> Connection:
    if username is None:
        username = os.getenv(ENV_WAREHOUSE_USERNAME)
    if password is None:
        password = os.getenv(ENV_WAREHOUSE_PASSWORD)
    if host is None:
        host = os.getenv(ENV_WAREHOUSE_HOST)
    if port is None and ENV_WAREHOUSE_PORT in os.environ:
        port = int(os.getenv(ENV_WAREHOUSE_PORT))
    if port is None:
        port = DEFAULT_WAREHOUSE_PORT
    if dbname is None:
        dbname = os.getenv(ENV_WAREHOUSE_DBNAME) or DEFAULT_WAREHOUSE_DBNAME

    if username is None or password is None or host is None:
        raise ValueError(
            "warehouse_username, warehouse_password, and warehouse_url must be specified"
        )

    if sslmode is None:
        sslmode = os.getenv(ENV_WAREHOUSE_SSLMODE) or DEFAULT_WAREHOUSE_SSLMODE
    if sslmode not in ("verify-ca", "verify-full"):
        LOG.warn(f"Accessing Benchling Warehouse using a discouraged SSL mode {sslmode}")

    return pg.connect(
        dbname=dbname,
        user=username,
        password=password,
        host=host,
        port=port,
        sslmode=sslmode,
    )


def select_one(conn: Connection, query: str, **kwargs: Any) -> Dict[str, Any]:
    """
    Executes a `query`, substituting any variables with `kwargs`, and returns a single row as
    a dict with keys corresponding to the column names.
    """
    cur = conn.cursor()
    cur.execute(query, kwargs)
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"No entry found for query '{query}' with values {kwargs}")
    return dict(zip((column[0] for column in cur.description), row))


def select_all(conn: Connection, query: str, **kwargs: Any) -> pd.DataFrame:
    """
    Executes a `query`, substituting any variables with `kwargs`, and returns a DataFrame
    containing the results.
    """
    cur = conn.cursor()
    cur.execute(query, kwargs)
    rows = cur.fetchall()
    columns = [x[0] for x in cur.description]
    return pd.DataFrame(rows, columns=columns)


def _add_genome(
    row: pd.Series, genomes: Dict[str, str], default_genome: Optional[str]
) -> pd.Series:
    species = row["species"]
    if species in genomes:
        row["genome_build"] = genomes[species]
    elif default_genome is not None:
        row["genome_build"] = default_genome
    else:
        raise ValueError("No genome found for species {species}")
    return row


def create_metasheet(
    conn: Connection, ctb_id: str, genomes: Dict[str, str], default_genome: Optional[str]
) -> pd.DataFrame:
    sample_df = select_all(conn, ASSAY_QUERY, ctb_id=ctb_id)
    sample_df = sample_df.apply(
        _add_genome, genomes=genomes, default_genome=default_genome, axis=1
    )
    tracking_dict = select_one(conn, TRACKING_QUERY, ctb_id=ctb_id)
    if len(sample_df) != tracking_dict["number_of_samples"]:
        raise ValueError(
            f"Number of samples in assay table ({len(sample_df)}) does not match "
            f"number of samples in tracking table ({tracking_dict['number_of_samples']})"
        )
    return sample_df


def create_metasheet_from_benchling(
    *,
    ctb_id: str,
    warehouse_username: Optional[str] = None,
    warehouse_password: Optional[str] = None,
    warehouse_host: Optional[str] = None,
    warehouse_sslmode: Optional[str] = None,
    genomes_json: Optional[Path] = None,
    default_genome: Optional[str] = None,
    output_file: Optional[Path] = None,
) -> None:
    """Creates a metasheet from Benchling metadata for a given CTB ID.

    Args:
        warehouse_username: The username to use when connecting to the Benchling warehouse. If
            not specified, then the environment variable `WAREHOUSE_USERNAME` will be used.
        warehouse_password: The password to use when connecting to the Benchling warehouse. If not
            specified, then the environment variable `WAREHOUSE_PASSWORD` will be used.
        warehouse_host: The hostname to use when connecting to the Benchling warehouse. If not
            specified, then the environment variable `WAREHOUSE_HOST` will be used if it is set,
            otherwise the default hostname.
        warehouse_sslmode: The SSL mode to use when connecting to the Benchling warehouse. If not
            specified, then the environment variable `WAREHOUSE_SSLMODE` will be used if it is set,
            otherwise the default SSL mode.
        genomes_json: A JSON file containing a mapping of species to genome build.
        default_genome: The default genome build to use if no genome is found for a given species.
        output_file: The path to the metasheet file to write.
    """
    genomes: Dict[str, str] = {}
    if genomes_json is not None:
        with open(genomes_json) as fd:
            genomes = json.load(fd)

    conn = warehouse_connect(
        username=warehouse_username,
        password=warehouse_password,
        host=warehouse_host,
        sslmode=warehouse_sslmode,
    )

    sample_df = create_metasheet(
        conn, ctb_id=ctb_id, genomes=genomes, default_genome=default_genome
    )
    sample_df.to_csv(output_file or sys.stdout, sep="\t", index=False)

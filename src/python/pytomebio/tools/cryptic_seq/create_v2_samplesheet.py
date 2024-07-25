import json
import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd


def remove_Ns(v1_index2: str, umi_cycles: int) -> str:
    """
    Remove trailing 'N' characters from the string v1_index2.
    Ensure the resulting string has a length equal to umi_cycles.
    """
    assert (
        v1_index2.count("N") == umi_cycles
    ), f"UMI length in index 1 ({v1_index2.count('N')}) != umi_cycles ({umi_cycles})"
    v2_index2 = v1_index2.rstrip("N")
    return v2_index2


def populate_cloud_data(v1_data: pd.DataFrame, ProjectName: str, umi_cycles: int) -> list:
    """
    Populate cloud data from v1 data.

    Parameters:
        v1_data (pandas.DataFrame): DataFrame containing v1 data.
        ProjectName (str): Name of the project.
        umi_cycles (int): Number of cycles for index 1.

    Returns:
        list of dict: List of dictionaries representing cloud data.

    This function populates cloud data from the provided v1 data. It first copies
    the input DataFrame `v1_data`. Then, it processes the 'index2' column to remove
    trailing 'N' characters using the remove_Ns function. Next, it adds the 'ProjectName'
    column with the specified project name. It creates the 'LibraryName' column by
    concatenating 'Sample_ID', 'index', and processed 'index2' columns with underscores.
    The column names are then converted to lowercase and renamed to match the expected
    format. Finally, the data is converted to a list of dictionaries with each dictionary
    representing a row of cloud data, and returned.
    """
    cloud_data = v1_data.copy()
    cloud_data["index2"] = cloud_data["index2"].apply(remove_Ns, umi_cycles=umi_cycles)
    cloud_data["ProjectName"] = ProjectName
    cloud_data["LibraryName"] = (
        cloud_data["Sample_ID"] + "_" + cloud_data["index"] + "_" + cloud_data["index2"]
    )

    cols = {
        "Sample_ID": "sample_id",
        "ProjectName": "project_name",
        "LibraryName": "library_name",
    }
    cloud_data = cloud_data.rename(columns=cols)[list(cols.values())]
    return cloud_data.to_dict(orient="records")


def populate_bclconvert_data(
    v1_data: pd.DataFrame,
    umi_cycles: int,
) -> list:
    """
    Populate bclconvert data from v1 data.

    Parameters:
        v1_data (pandas.DataFrame): DataFrame containing v1 data.
        umi_cycles (int): Number of cycles for index 1.

    Returns:
        list of dict: List of dictionaries representing bclconvert data.

    This function populates bclconvert data from the provided v1 data. It extracts
    columns 'Sample_ID', 'index', and 'index2' from v1_data. The 'index2' column
    values are processed to remove trailing 'N' characters using the remove_Ns
    function. The column names are converted to lowercase. Finally, the data is
    converted to a list of dictionaries with each dictionary representing a row
    of bclconvert data, and returned.
    """
    cols = ["Sample_ID", "index", "index2"]
    bclconvert_data = v1_data[cols].copy()
    bclconvert_data["index2"] = bclconvert_data["index2"].apply(remove_Ns, umi_cycles=umi_cycles)
    bclconvert_data.columns = bclconvert_data.columns.str.lower()
    return bclconvert_data.to_dict(orient="records")


def get_umi_cycles(v2_metasheet: dict) -> int:
    """
    Get the number of cycles for the UMI (Unique Molecular Identifier) in the v2 metasheet.

    Parameters:
        v2_metasheet (dict): Dictionary representing the v2 metasheet.

    Returns:
        int: Number of cycles for the UMI.

    This function retrieves the number of cycles for the Unique Molecular Identifier (UMI)
    from the v2 metasheet. It expects a dictionary representing the v2 metasheet as input.
    The UMI cycles information is extracted from the 'override_cycles' field in the 'reads'
    section of the v2 metasheet. The UMI cycles are assumed to be specified as the third
    component separated by semicolons in the 'override_cycles' string. If the UMI cycles
    are prefixed with 'U', it extracts the numeric part after 'U' and returns it as an integer.
    """
    override_cycles = v2_metasheet["bclconvert_settings"]["override_cycles"]
    index_2_cycles = override_cycles.split(";")[2]
    return int(index_2_cycles.split("U")[-1])


def v2_metasheet_template(v2_template: Path) -> dict:
    """
    Load v2 metasheet template from JSON file.

    Parameters:
        v2_template (str): Path to the v2 metasheet template JSON file.

    Returns:
        dict: Dictionary representing the v2 metasheet template.
    """
    # Open and load the JSON file
    with open(v2_template, "r") as file:
        return json.load(file)


def v1_samplesheet_data(v1_samplesheet: Path) -> pd.DataFrame:
    """
    Read from a v1 Illumina samplesheet and return the [Data] section.

    Parameters:
        v1_samplesheet: Path to the v1 Illumina samplesheet file.

    Returns:
        pd.DataFrame: DataFrame containing the data from the samplesheet.

    """
    with open(v1_samplesheet, "r") as file:
        in_data_section = False
        data_lines = []

        for line in file:
            if line.strip() == "[Data]":
                in_data_section = True
                continue
            elif line.strip().startswith("[") and line.strip() != "[Data]":
                if in_data_section:
                    break

            if in_data_section:
                data_lines.append(line.strip())

        if not data_lines:
            raise ValueError("No [Data] section found or it is empty")

        # Combine the collected lines into a single string
        data_str = "\n".join(data_lines)

        # Use StringIO to simulate a file object for pandas
        data_io = StringIO(data_str)

        # Read the data into a pandas DataFrame
        return pd.read_csv(data_io)


def create_v2_samplesheet(
    *,
    v1_samplesheet: Path,
    v2_template: Path,
    output_file: Path,
    ctb_id: str,
) -> None:
    """
    Creates the Illumina sample sheet in v2 format to demultiplex
    using bcl-convert.

    Args:
        v1_samplesheet (Path): Path to the v1 Illumina samplesheet file.
        v2_template (Path): Path to the v2 metasheet template JSON file.
        ctb_id (str): CTB ID used for the run name and project name.
        output_file (Path): Path to the output v2 samplesheet file.

    This function generates a v2 samplesheet for use with the bcl-convert tool.
    It takes a v1 Illumina samplesheet, a v2 metasheet template, and a CTB ID as inputs.
    The v1 samplesheet is read and processed to extract relevant information.
    The v2 metasheet template is loaded and modified with the CTB ID as the run name.
    The bclconvert and cloud data sections of the v2 metasheet are populated using
    information from the v1 samplesheet.
    The generated v2 samplesheet is written to the specified output file.
    """

    v1_data = v1_samplesheet_data(v1_samplesheet)
    v2_metasheet = v2_metasheet_template(v2_template)
    umi_cycles = get_umi_cycles(v2_metasheet)
    v2_metasheet["header"]["run_name"] = ctb_id
    v2_metasheet["bclconvert_data"] = populate_bclconvert_data(v1_data, umi_cycles=umi_cycles)
    v2_metasheet["cloud_data"] = populate_cloud_data(
        v1_data, ProjectName=ctb_id, umi_cycles=umi_cycles
    )

    with open("v2_tmp.json", "w") as outfile:
        json.dump(v2_metasheet, outfile)
    command = ["v2-samplesheet-maker", "v2_tmp.json", str(output_file)]
    subprocess.check_call(command)

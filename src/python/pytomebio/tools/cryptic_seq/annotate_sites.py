import csv
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List
from typing import Optional

import pandas as pd
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pytomebio import core
from pytomebio.pipeline.samples import FileGrouping
from pytomebio.pipeline.samples import SampleWithSites


def fetch_sequence(
    sample: SampleWithSites, attP_half_site_len: int, attB_half_site_len: int
) -> pd.DataFrame:
    """
    Takes a sites.txt, finds the boundaries of the full cryptic site
    and fetches the sequence from the fasta file.

    Returns a dataframe with the following added columns to the original ones in sites.txt:
        seq_start
        seq_end
        strand
        seq
        genome_dinucleotide
    """

    df = pd.read_csv(sample.metrics_file, sep="\t")

    # Find boundaries of full cryptic att based on attachment_site
    # TO DO: check if seq_start < 1, if it is fix before bedtools call and add Ns to the sequence
    df["seq_start"] = df.apply(
        lambda x: x["position"] - attP_half_site_len
        if "attP" in x["attachment_site"]
        else x["position"] - attB_half_site_len,
        axis=1,
    )
    df["seq_end"] = df.apply(
        lambda x: x["position"] + 2 + attP_half_site_len
        if "attP" in x["attachment_site"]
        else x["position"] + 2 + attB_half_site_len,
        axis=1,
    )
    df["strand"] = df["positive_strand"].apply(lambda x: "+" if x else "-")
    df["name"] = (
        df["attachment_site"]
        + ":"
        + df["reference_name"]
        + ":"
        + df["seq_start"].astype(str)
        + ":"
        + df["seq_end"].astype(str)
        + ":"
        + df["strand"].astype(str)
    )

    # make sure name is unique
    assert not any(df["name"].duplicated()), f"Duplicated sites in {sample.metrics_file}"

    # fetch sequence using bedtools getfasta with options
    # s: force strandness, tab: tab-delimited format, name:
    cols = ["reference_name", "seq_start", "seq_end", "name", "count", "strand"]
    bed = pybedtools.BedTool.from_dataframe(df[cols])
    seqs = bed.sequence(fi=str(sample.ref_fasta), name=True, s=True, tab=True)
    seqs = pd.read_csv(seqs.seqfn, sep="\t", names=["name", "seq"])
    seqs["seq"] = seqs["seq"].str.upper()
    seqs["name"] = seqs["name"].str.split("::").str[0]
    assert len(seqs) == len(df)

    # merge sites df with seqs to add column seq
    df = pd.merge(df, seqs).drop(columns=["name"])
    # find the genome dinucleotide from seq
    df["genome_dinucleotide"] = df.apply(
        lambda x: x.seq[attP_half_site_len : attP_half_site_len + 2]
        if "attP" in x["attachment_site"]
        else x.seq[attB_half_site_len : attB_half_site_len + 2],
        axis=1,
    )
    return df


def recombination_type(row: pd.Series) -> str:
    donor_dn = row["attachment_site"][-2:].upper()
    assert all(x in "ACGT" for x in donor_dn)
    if donor_dn == row["genome_dinucleotide"]:
        return "canonical"
    else:
        return "crosstalk"


parse_group = core.key_value_parser(r'(\w+)\s+"([^"]*)"\s*;')


def get_threat_tier(row: pd.Series) -> pd.Series:
    """
    Get the threat tier based on list of features in a row.
    If there are more than one overlapping feature (expected if gene contains multiple transcript
    isoforms), the following criteria is applied:
        * CDS as feauture in any of the transcript isoforms -> tier I
        * no CDS but any of the transcript isoform is protein coding-> tier II
        * overlaps non-protein coding gene only -> tier III
        * no overlapping gene -> tier IV
    """
    assert len(row.feature) == len(row.group)

    has_overlap = False
    has_gene = False
    for feature, group in zip(row.feature, row.group):
        if feature == "." and group == ".":
            continue
        else:
            has_overlap = True

        if feature == "CDS":
            return "I"
        elif feature == "gene":
            has_gene = True
            annotations = parse_group(group)
            if annotations["gene_biotype"] == "protein_coding":
                return "II"

    if has_gene:
        return "III"
    elif not has_overlap:
        return "IV"
    else:
        raise Exception("Row does not contain a CDS or gene feature: {row}")


def overlapping_gene(row: pd.Series) -> Optional[str]:
    """
    Returns information about the genes that overlap with each cryptic site
    in the format of [<gene_id>,<strand>,<gene_biotype>;]
    """
    genes = [
        (parse_group(group), strand)
        for feature, group, strand in zip(row.feature, row.group, row._strand)
        if feature == "gene"
    ]

    if len(genes) == 0:
        return None

    return "; ".join(
        [f"{group['gene_id']},{strand},{group['gene_biotype']}" for group, strand in genes]
    )


def overlapping_feature(row: pd.Series) -> Optional[str]:
    """
    Returns the specific gene feature, if any, that overlaps with a cryptic site.

    Considers all features (from multiple transcript isoforms if present in the gene set)
    and returns the corresponding feature with the following criteria:
        'CDS' among features -> 'CDS'
        no 'CDS' but 'exon' is among features -> 'intron'
        else if biotype is protein coding (mRNA) -> 'UTR' (untranslated region)
        else if 'exon' is among features (non-coding gene) -> 'non_coding_exon'

    TODO: double check other biotypes would make it to the last category
    """
    has_gene = False
    has_exon = False
    has_utr = False
    has_cds = False

    for feature, group in zip(row.feature, row.group):
        if feature == "gene":
            has_gene = True
        elif feature == "CDS":
            has_cds = True
        elif feature == "exon":
            has_exon = True
            annotations = parse_group(group)
            if annotations["transcript_biotype"] == "mRNA":
                has_utr = True

    if not has_gene:
        # no gene overlap
        return None
    elif not has_exon:
        # gene overlap but no exon, must be intron
        return "intron"
    elif has_cds:
        # CDS overlap
        return "CDS"
    elif has_utr:
        # no CDS, exon overlap in protein_coding must be UTR
        return "UTR"
    else:
        # exon in transcript_biotype other than mRNA
        return "non_coding_exon"


def same_strand(row: pd.Series) -> str:
    """
    Check whether the cryptic site and gene are on the same strand.
    Returns 'both' if a site overlaps with 2 genes genes, each on a different strand.
    """
    if not row.overlapping_gene:
        return None

    strands = set(gene.split(",")[1] for gene in row.overlapping_gene.split("; "))

    if "-" in strands and "+" in strands:
        return "both"
    elif row.strand in strands:
        return "yes"
    elif len(strands) > 0:
        return "no"
    else:
        raise Exception(f"Cannot determine strand from: {overlapping_gene}")


def cosmic_genes(df: pd.DataFrame, gtf: Path) -> pd.DataFrame:
    """
    Annotates the input dataframe with the following columns from the COSMIC database:
    """
    raise NotImplementedError()


def gene_overlaps(df: pd.DataFrame, gtf: Path) -> pd.DataFrame:
    """
    Computes the overlaps between the cryptic sites and the gene annotation provided as GTF.

    Returns the input dataframe with the following added columns:
        threat_tier: 4 tier category for cryptic sites based on the overlapping gene feature
            Tier I: coding regions only
            Tier II: non-coding regions (exonic UTR or intronic) of coding genes
            Tier III: exonic or intronic regions of non-coding genes
            Tier IV: all regions not in Tier I-III
        overlapping_gene: information about overlapping genes in the format of
            <gene_id>,<strand>,<gene_biotype>;
        overlapping_feature: specific gene feature overlapped by the cryptic site
        same_strand: whether the cryptic site and the gene are in the same strand,
            coulb be "both" if the site overlaps with gene on both strands

    Runs pybedtools.intersect
    """

    # create name based on genomic coordinates, it will be used below to merge with annotated df
    df["name"] = (
        df["reference_name"]
        + ":"
        + df["seq_start"].astype(str)
        + ":"
        + df["seq_end"].astype(str)
        + ":"
        + df["strand"]
    )
    original_df = df.copy()

    # keep columns for bed only and remove redundancy in sites
    # sites can be duplicated if multiple donors were used
    df["score"] = 0
    cols_for_bed = ["reference_name", "seq_start", "seq_end", "name", "score", "strand"]
    df = df[cols_for_bed].drop_duplicates()

    # make sure name is unique
    duplicated_rows = df[df["name"].duplicated()]
    assert duplicated_rows.empty, f"Duplicated sites:\n{duplicated_rows}"

    # Convert numeric columns to integers
    numeric_cols = df.select_dtypes(include=["number"]).columns
    df[numeric_cols] = (
        df[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(df[numeric_cols]).astype(int)
    )

    # sites as bedtool object
    a_bed = pybedtools.BedTool.from_dataframe(df[cols_for_bed])
    b_gtf = pybedtools.BedTool(gtf)
    # overlap sites and gene set
    # run bedtools intersect with options
    # loj: left outer join; wb: write original entry in b (gtf); sorted: a (sites) is sorted
    names = cols_for_bed + [
        "seqname",
        "source",
        "feature",
        "_start",
        "_end",
        "_score",
        "_strand",
        "frame",
        "group",
    ]
    intersected_df = (
        a_bed.sort()
        .intersect(b=b_gtf, sorted=True, wb=True, loj=True, nonamecheck=True)
        .to_dataframe(names=names)
    )

    # group by coordinates and collapse columns feature and group as lists
    # each row now corresponds to a unique site, and contains the list of overlapping features
    by_fields = ["reference_name", "seq_start", "seq_end", "strand", "name"]
    grouped_data = (
        intersected_df.groupby(by_fields)
        .agg({"feature": list, "group": list, "_strand": list})
        .reset_index()
    )
    # make sure all sites are here
    assert len(a_bed) == len(grouped_data)

    # get thread tier for each site
    # the tier is based on the overlaps returned by BedTool.intersect. Briefly,
    # if a site overlaps with CDS it returns Tier I
    # else if overlaps with a coding gene, returns Tier II
    # else if overlaps with a gene whose biotype is other than protein_coding, returns Tier III
    # else if site overlaps no features in the gtf file, returns IV
    grouped_data["threat_tier"] = grouped_data.apply(get_threat_tier, axis=1)
    # get overlapping gene and the corresponding feature (CDS,intron,UTR,non_coding_exon)
    grouped_data["overlapping_gene"] = grouped_data.apply(overlapping_gene, axis=1)
    grouped_data["overlapping_feature"] = grouped_data.apply(overlapping_feature, axis=1)
    grouped_data["same_strand"] = grouped_data.apply(same_strand, axis=1)

    # drop temp colums
    grouped_data = grouped_data.drop(columns=["feature", "group", "_strand"])

    # merge with original df
    df = pd.merge(original_df, grouped_data)
    return df


def collate_sites(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collate sites by position and strand.

    Returns a dataframe with unique list of cryptic sites:
        sample_donor_count: keeps track of individual recombination events in the format of
            <sample>:<donor>:<count>
        nsamples: number of samples where a cryptic site has ben detected
            (does not consider multiple event in the same sample ie. multiple donors)
    """

    # aggregate sites by genomic coordinate
    unique_sites_group_cols = ["reference_name", "position", "seq_start", "seq_end", "strand"]
    unique_sites_df = (
        df.groupby(unique_sites_group_cols)
        .agg(
            {
                "sample": list,
                "attachment_site": list,
                "count": list,
                # these should be 1 element sets if the same fasta was used for all samples
                "seq": set,
                "genome_dinucleotide": set,
                # these should be 1 element sets if the same gtf was used
                "threat_tier": set,
                "overlapping_gene": set,
                "overlapping_feature": set,
                "same_strand": set,
            }
        )
        .reset_index()
    )
    # List of columns that are of type 'set'
    set_columns = [
        "seq",
        "genome_dinucleotide",
        "threat_tier",
        "overlapping_gene",
        "overlapping_feature",
        "same_strand",
    ]

    # Iterate over each set column and join the elements with commas
    for col in set_columns:
        unique_sites_df[col] = unique_sites_df[col].apply(
            lambda x: ",".join(x) if x != {None} else None
        )

    # count number of samples
    unique_sites_df["nsamples"] = unique_sites_df["sample"].apply(lambda x: len(set(x)))

    # keep track of each recombination event observed for a particular site
    # the new column sample_donor_count will be a list of events
    # with information about sample, donor, and count
    def custom_combine(row: pd.Series) -> pd.Series:
        zipped = zip(row["sample"], row["attachment_site"], row["count"])
        return ";".join(f"{sample}:{donor}:{count}" for sample, donor, count in zipped)

    unique_sites_df["sample_donor_count"] = unique_sites_df.apply(custom_combine, axis=1)
    unique_sites_df = unique_sites_df.drop(columns=["sample", "attachment_site", "count"])

    # compute mean counts across all events
    all_events_count_df = (
        df.groupby(unique_sites_group_cols)
        .agg(
            {
                "count": "mean",
            }
        )
        .reset_index()
    )

    # compute mean counts across events of the same recombination type: canonical and cross-talk
    recombination_group_cols = [
        "reference_name",
        "position",
        "seq_start",
        "seq_end",
        "strand",
        "recombination",
    ]
    recombination_count_df = (
        df.groupby(recombination_group_cols)
        .agg(
            {
                "count": "mean",
            }
        )
        .reset_index()
    )

    # Pivot the DataFrame
    recombination_count_df = recombination_count_df.pivot_table(
        index=unique_sites_group_cols, columns="recombination", values="count"
    ).reset_index()

    # merge all events and recombination dfs
    mean_count_df = (
        pd.merge(all_events_count_df, recombination_count_df)
        .rename(
            columns={
                "count": "mean_count",
                "canonical": "mean_count_canonical",
                "crosstalk": "mean_count_crosstalk",
            }
        )
        .reset_index(drop=True)
    )

    # merge mean counts with unique sites dfs
    df = pd.merge(mean_count_df, unique_sites_df).reset_index(drop=True)

    # transform start coordinates to 1-based
    df["dinucleotide_position"] = df["position"] + 1

    # final columns selection and order
    column_order = [
        "reference_name",
        "dinucleotide_position",
        "strand",
        "genome_dinucleotide",
        "mean_count",
        "mean_count_canonical",
        "mean_count_crosstalk",
        "seq_start",
        "seq_end",
        "seq",
        "threat_tier",
        "overlapping_gene",
        "overlapping_feature",
        "same_strand",
        "nsamples",
        "sample_donor_count",
    ]

    return df[column_order]


def parse_fastmap_results(path: Path) -> pd.DataFrame:
    with open(path) as results:
        reader = csv.reader(results, delimiter="\t")
        name = None
        seq_len = None
        num_hits = 0
        rows = []
        for line in reader:
            assert line[0] == "SQ"
            name = line[1]
            seq_len = int(line[2])

            line = next(reader)
            while line[0] != "//":
                assert line[0] == "EM"
                start = int(line[1])
                end = int(line[2])
                if end - start == seq_len:
                    num_hits += int(line[3])
                line = next(reader)
                # there are sometimes blank lines
                while len(line) == 0:
                    line = next(reader)

            rows.append((name, num_hits))
            name = None
            seq_len = None
            num_hits = 0

    return pd.DataFrame(rows, columns=["name", "genome_hits"])


def fastmap(query_fasta: Path, bwa_index_base: Path) -> pd.DataFrame:
    output_file = NamedTemporaryFile(mode="w", delete=False)
    output_path = Path(output_file.name)
    command = ["bwa", "fastmap", str(bwa_index_base), str(query_fasta)]
    try:
        subprocess.check_call(command, stdout=output_file)
        output_file.close()
        return parse_fastmap_results(output_path)
    finally:
        try:
            # output_path.unlink()
            print(output_path)
        except IOError:
            pass


def find_genome_exact_matches(df: pd.DataFrame, bwa_index_base: Path) -> pd.DataFrame:
    df["name"] = (
        df["reference_name"]
        + ":"
        + df["seq_start"].astype(str)
        + ":"
        + df["seq_end"].astype(str)
        + ":"
        + df["strand"]
    )
    query_fasta_file = NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    query_fasta_path = Path(query_fasta_file.name)
    try:
        # write query fasta file
        num_written = SeqIO.write(
            (
                SeqRecord(Seq(row["seq"]), id=row["name"], description="")
                for _, row in df.iterrows()
            ),
            query_fasta_file,
            "fasta",
        )
        query_fasta_file.close()
        assert num_written == len(df)
        genome_hits = fastmap(query_fasta_path, bwa_index_base)
        assert df.shape[0] == genome_hits.shape[0], f"{df.shape[0]} != {genome_hits.shape[0]}"
        return pd.merge(df, genome_hits, how="inner").drop(columns=["name"])
    finally:
        try:
            query_fasta_path.unlink()
        except IOError:
            pass


def explode_by_gene(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a dataframe containing each gene overlapping an integration site.
    If an integration site overlaps more than one gene, each will be present in different rows.
    """
    # skip sites that do not overlap with genes
    df = df[df["overlapping_gene"].notna()].copy()

    df["overlapping_gene"] = df["overlapping_gene"].str.split(";")
    df = df.explode("overlapping_gene")
    df[["gene_name", "gene_strand", "gene_biotype"]] = df["overlapping_gene"].str.split(
        ",", expand=True
    )
    cols = [
        "gene_name",
        "reference_name",
        "dinucleotide_position",
        "strand",
        "gene_strand",
        "same_strand",
        "gene_biotype",
        "genome_dinucleotide",
        "threat_tier",
        "overlapping_feature",
        "nsamples",
        "mean_count_canonical",
        "mean_count_crosstalk",
    ]
    df = df[cols].sort_values(by="gene_name", ascending=True)
    return df


def annotate_sites(
    *,
    in_yml: Optional[Path] = None,
    in_json: Optional[Path] = None,
    ref_dir: Optional[Path] = None,
    match_index_base: Optional[Path] = None,
    output: Path = Path("sites.annotated.xlsx"),
    file_grouping: FileGrouping = FileGrouping.group_name,
    attP_half_site_len: int = 22,
    attB_half_site_len: int = 25,
) -> None:
    """Annotates sites and compiles output tables into a multi-sheet excel file.

    The input YAML file must be the same config file used in the CHANGE-Seq and Cryptic-Seq
    Snakemake pipelines. Please refer to the top-level README for the current format. Alternately,
    an input JSON file with a list of sample objects may be provided.

    Args:
        in_yml: the input configuration YAML used in Cryptic-Seq Snakemake
                pipeline.
        in_json: the input JSON file with a list of sample objects.
        ref_dir: the path to the directory with reference FASTA files, if the reference paths in
            the input file are relative.
        match_index_base: the path of the BWA index base to use when finding matches for site
            sequences in the genome.
        output: output excel file, contains annotated output tables.
        file_grouping: how files are organized in the working directory.
        attP_half_site_len: length of the genomic half attachment site (without dinucleotide)
                detected with an attP donor, that is a cryptic attB site.
        attB_half_site_len: length of the genomic half attachment site (without dinucleotide)
                detected with an attB donor, that is a cryptic attP site.
    """
    # Load in the samples using the pipeline config
    samples: List[SampleWithSites] = None
    if in_yml:
        samples = SampleWithSites.from_yml(path=in_yml, grouping=file_grouping, ref_dir=ref_dir)
    elif in_json:
        samples = SampleWithSites.from_json(path=in_json, grouping=file_grouping, ref_dir=ref_dir)
    else:
        raise ValueError("in_yml or in_json must be provided")

    # Create an Excel writer object
    writer = pd.ExcelWriter(output, engine="openpyxl")

    # Loop samples, for each site identified
    # fetch sequence, find genome dinucleotide and define recombination type
    list_of_df = []
    for sample in samples:
        # add full cryptic site sequence and genome dinucleotide
        df = fetch_sequence(sample, attP_half_site_len, attB_half_site_len)

        # recombination type based on match between donor and genome dinucleotide
        df["recombination"] = df.apply(recombination_type, axis=1)

        # compute gene overlaps
        df = gene_overlaps(df, sample.gene_set_gtf)

        # transform start coordinates to 1-based
        df["dinucleotide_position"] = df["position"] + 1
        df["seq_start"] += 1

        df["group"] = sample.group
        df["sample"] = sample.name

        column_order = [
            "group",
            "sample",
            "reference_name",
            "dinucleotide_position",
            "strand",
            "attachment_site",
            "genome_dinucleotide",
            "recombination",
            "seq_start",
            "seq_end",
            "count",
            "seq",
            "threat_tier",
            "overlapping_gene",
            "overlapping_feature",
            "same_strand",
        ]

        # add sample to excel
        df[column_order].to_excel(writer, sheet_name=sample.name, index=False)
        # add sample to list of dfs
        list_of_df.append(df)

    # concat all sites
    df = pd.concat(list_of_df)
    # collate sites by chr, position, and strand
    # keeping track of observed recombination events and compute mean counts
    df = collate_sites(df)
    # find exact matches in the genome for each cryptic site
    if match_index_base:
        df = find_genome_exact_matches(df, match_index_base)
    # TODO: overlap with COSMIC genes
    df.to_excel(writer, sheet_name="sites_collated", index=False)

    # aggregate by gene
    by_gene_df = explode_by_gene(df)
    by_gene_df.to_excel(writer, sheet_name="genes", index=False)

    # Close the Excel file
    writer.close()

    pybedtools.helpers.cleanup(verbose=False, remove_all=False)

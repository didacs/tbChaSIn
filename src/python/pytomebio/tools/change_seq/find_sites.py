from attr import frozen
from collections import Counter
from fgpyo.util.metric import Metric
from pathlib import Path
from samwell import sam
from samwell.itertools import peekable
from typing import Dict
from typing import Optional
from typing import List


@frozen
class FindSitesMetric(Metric["FindSitesMetric"]):
    reference_name: str
    position: int
    count: int
    site_to_count: Optional[Dict[int, int]]


def aggregate_sites(sites: Dict[str, Counter], slop: int) -> List[FindSitesMetric]:
    """Aggregates sites within the slop away from each other.

    The position will be set as the site with the highest read count.

    Args:
        sites: the sites to aggregate
        slop: the maximum distance for two sites to be aggregated
    """
    # the list of aggregated sites
    agg_sites: List[FindSitesMetric] = []
    for reference_name, counter in sites.items():

        # running stats about a given aggregated site
        prev_site_pos: int = -abs(slop) - 1  # the position of the most recent site encountered
        running_count: int = 0  # the running count of read pairs for the aggregated site
        cur_sites: Dict[int, int] = {}  # the positions and counts for sites being aggregated
        max_count_pos: int = -1  # the position of the site with the highest count in the set
        max_count: int = -1  # the count of read pairs for the site with the highest count
        for pos, count in sorted(counter.items(), key=lambda tup: tup[0]):

            if pos - prev_site_pos > slop:  # new site
                if prev_site_pos > 0:  # add the previous site
                    site = FindSitesMetric(
                        reference_name=reference_name,
                        position=max_count_pos,
                        count=running_count,
                        site_to_count=cur_sites,
                    )
                    agg_sites.append(site)
                # re-initialize
                max_count = -1
                cur_sites = {}

            # add the current position to the site being aggregated
            prev_site_pos = pos
            running_count += count
            cur_sites[pos] = count
            if max_count < count:
                max_count = count
                max_count_pos = pos

        # add the last site found
        assert prev_site_pos > 0
        site = FindSitesMetric(
            reference_name=reference_name,
            position=max_count_pos,
            count=running_count,
            site_to_count=cur_sites,
        )
        agg_sites.append(site)

    return agg_sites


def find_sites(
    *,
    in_bam: Path,
    out_txt: Path,
    intra_read_slop: int = 5,
    inter_site_slop: int = 5,
    min_mapq: int = 20
) -> None:
    """Finds integration sites by examining read pair alignments.

    The read pairs are assumed to have the leading attachment site trimmed, without the overhang
    trimmed (see tomebio-tools tomechange-seq fastq-to-bam).

    The integration site should occur at the first sequenced base in each read in a read pair. Thus
    each read pair is examined and a site is called if the genomic position of the first sequenced
    base in each read pair is within the given intra-read-slop.  Next, the counts for each
    integration site in the genome is aggregated.  Finally, sites within inter-site-slop distance
    from each other are aggregated.

    For sites that are aggregated, the position is the site with the highest read count.

    Args:
        in_bam: path to the input BAM either sorted by queryname or grouped by query.
        out_txt: the list of sites found
        intra_read_slop: the maximum absolute difference between the first base sequenced within
                         a read pair.  Typically, at least 2bp to account for the overhang in the
                         attachment site, which is not trimmed.
        inter_site_slop: aggregate sites that are within this distance of each other.  The counts
                         are summed, and the coordinate is the site with the largest count.  Set to
                         zero to not aggregate counts.
        min_mapq: the minimum mapping quality to consider for a read pair
    """

    with sam.reader(in_bam, file_type=sam.SamFileType.BAM) as reader:
        assert (
            reader.header["HD"].get("SO") == "queryname"
            or reader.header["HD"].get("GO") == "query"
        ), "Input BAM must be queryname sorted or query grouped"

        # Stores the raw sites, unaggregated.  Mapping from contig name to a Counter, where each
        # Counter stores the count of reads for a given site's position.
        raw_sites: Dict[str, Counter] = {}

        # Read in all records with the same name
        piter = peekable(iter(reader))
        while piter.can_peek():
            first = next(piter)
            records = [first]
            while piter.can_peek() and piter.peek().query_name == first.query_name:
                records.append(next(piter))

            if len(records) <= 1:
                for record in records:
                    print(record)
            read1 = next(
                r for r in records if r.is_read1 and not r.is_secondary and not r.is_supplementary
            )
            read2 = next(
                r for r in records if r.is_read2 and not r.is_secondary and not r.is_supplementary
            )

            if read1.is_unmapped or read2.is_unmapped:
                continue
            if read1.reference_name != read2.reference_name:
                continue
            if read1.mapping_quality < min_mapq or read2.mapping_quality < min_mapq:
                continue

            # get the position of the first base _sequenced_
            read1_pos = read1.reference_start if read1.is_forward else read1.reference_end
            read2_pos = read1.reference_start if read1.is_forward else read1.reference_end
            pos_diff = abs(read1_pos - read2_pos)
            if pos_diff > intra_read_slop:
                continue
            pos = (read1_pos + read2_pos) // 2

            if read1.reference_name not in raw_sites:
                raw_sites[read1.reference_name] = Counter()
            counter = raw_sites[read1.reference_name]
            counter[pos] = counter[pos] + 1

        # Aggregate sites
        all_sites: List[FindSitesMetric] = aggregate_sites(sites=raw_sites, slop=inter_site_slop)

        # Sort the sites by reference_name
        reference_name_to_int: Dict[str, int] = {
            name: i for i, name in enumerate(reader.header.references)
        }
        all_sites = sorted(
            all_sites, key=lambda m: (reference_name_to_int[m.reference_name], m.position)
        )

        # Write the metrics
        FindSitesMetric.write(out_txt, *all_sites)

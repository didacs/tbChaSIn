import enum
import itertools
from abc import ABC
from abc import abstractmethod
from collections import Counter
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple

from attr import define
from attr import frozen
from fgpyo import sam
from fgpyo.sam import Template
from fgpyo.util.metric import Metric


@define
class FindSitesMetric(Metric["FindSitesMetric"]):
    """Collates read counts for an integration site.

    Attributes:
        reference_name: the genomic contig or reference name
        position: the genomic position (1-based)
        attachment_site: the name of the attachment site (e.g. attP or attB)
        positive_strand: true if the left side of the attachment site integrates on the positive
                         strand of the genome, false otherwise
        left_to_count: mapping of genomic position to template count for observations where
                       just the integration of left side of the attachment site was observed in a
                       single template
        right_to_count: mapping of genomic position to template count for observations where
                       just the integration of right side of the attachment site was observed in a
                       single template
        both_to_count: mapping of genomic position to template count for observations where
                       the integration of both sides of the attachment site was observed in a
                       single template
        count: the total count of templates observing the integration of this attachment site (sum
               over left, right, and both)
    """

    reference_name: str
    position: int
    attachment_site: str
    positive_strand: bool
    left_to_count: Counter
    right_to_count: Counter
    both_to_count: Counter
    count: int


@enum.unique
class Side(enum.Enum):
    """The side of the attachment site observed in a single template (read pair) for an
    integration site.

    Change-seq observes both sides in the same template, while Cryptic-seq observes only one side.
    """

    Left = "left"
    Right = "right"
    Both = "both"


@frozen
class SiteKey:
    """Key for `FindSitesMetric`s that should to gather template counts in
    `SitesGenerator.generate`.

    Attributes:
        reference_name: the genomic contig or reference name
        position: the genomic position (1-based)
        attachment_site: the name of the attachment site (e.g. attP or attB)
        positive_strand: true if the left side of the attachment site integrates on the positive
                         strand of the genome, false otherwise
    """

    reference_name: str
    position: int
    attachment_site: str
    positive_strand: bool


def aggregate_sites(sites: List[FindSitesMetric], slop: int) -> List[FindSitesMetric]:
    """Aggregates sites within the slop away from each other.

    The position will be set as the site with the highest read count.

    Args:
        sites: the sites to aggregate
        slop: the maximum distance for two sites to be aggregated
    """
    # the list of aggregated sites
    agg_sites: List[FindSitesMetric] = []

    # sort and then group by reference name, attachment site, and positive strand
    def sort_key(site: FindSitesMetric) -> Tuple[str, str, bool]:
        return site.reference_name, site.attachment_site, site.positive_strand

    for _, group in itertools.groupby(sorted(sites, key=sort_key), key=sort_key):
        cur_site: Optional[FindSitesMetric] = None
        for site in sorted(group, key=lambda s: s.position):
            if cur_site is None or site.position - cur_site.position > slop:  # new site
                if cur_site is not None:
                    counter = (
                        cur_site.left_to_count + cur_site.right_to_count + cur_site.both_to_count
                    )
                    cur_site.position = counter.most_common(1)[0][0]
                    agg_sites.append(cur_site)
                # re-initialize
                cur_site = site
            else:
                cur_site.position = site.position
                cur_site.left_to_count.update(site.left_to_count)
                cur_site.right_to_count.update(site.right_to_count)
                cur_site.both_to_count.update(site.both_to_count)
                cur_site.count += site.count

        if cur_site is not None:
            counter = cur_site.left_to_count + cur_site.right_to_count + cur_site.both_to_count
            cur_site.position = counter.most_common(1)[0][0]
            agg_sites.append(cur_site)

    return agg_sites


class SitesGenerator(ABC):
    @abstractmethod
    def get_key(self, template: Template) -> Optional[SiteKey]:
        """Returns the key used to identify an integration site, None if not to call an integration
        site.

        The reads are assumed to be mapped to the same contig.

        Args:
            template: the template (read pair)
        """
        pass

    @abstractmethod
    def get_side(self, template: Template) -> Optional[Side]:
        """Returns the "side" of the integration that was observed.

        Args:
            template: the template (read pair)
        """
        pass

    def generate(
        self,
        in_bam: Path,
        inter_site_slop: int,
        min_mapq: int,
        use_duplicates: bool,
        same_reference: bool,
    ) -> List[FindSitesMetric]:
        """Generates a list of integration sites by examining read pairs in the given BAM file.

        Args:
            in_bam: the path to the input BAM file
            inter_site_slop: aggregate sites that are within this distance of each other.  The
                             counts are summed, and the coordinate is the site with the largest
                             count.  Set to zero to not aggregate counts.
            min_mapq: the minimum mapping quality to consider for a read pair
            use_duplicates: true to ignore the duplicate flag on read pairs, false to not count
                            duplicate reads.
            same_reference: true to require hte primary alignments to be mapped to the same contig
        """
        with sam.reader(in_bam, file_type=sam.SamFileType.BAM) as reader:
            assert (
                reader.header["HD"].get("SO") == "queryname"  # type: ignore
                or reader.header["HD"].get("GO") == "query"  # type: ignore
            ), "Input BAM must be queryname sorted or query grouped"

            # First, we want to collate sites (and their counts) that have the following the same:
            # - reference name
            # - genomic position
            # - attachment site name (e.g. attP or attB)
            # - integration strand
            raw_sites: Dict[SiteKey, FindSitesMetric] = {}

            template_iterator: Iterator[Template] = Template.iterator(reader)
            for template in template_iterator:
                read1 = template.r1
                read2 = template.r2

                if read1.is_unmapped or read2.is_unmapped:
                    continue
                if same_reference and read1.reference_name != read2.reference_name:
                    continue
                if read1.mapping_quality < min_mapq or read2.mapping_quality < min_mapq:
                    continue
                if not use_duplicates and (read1.is_duplicate or read2.is_duplicate):
                    continue

                # get the position of the first base _sequenced_
                key = self.get_key(template)
                if key is None:
                    continue

                # add the key if not present
                if key not in raw_sites:
                    raw_sites[key] = FindSitesMetric(
                        reference_name=key.reference_name,
                        position=key.position,
                        attachment_site=key.attachment_site,
                        positive_strand=key.positive_strand,
                        left_to_count=Counter(),
                        right_to_count=Counter(),
                        both_to_count=Counter(),
                        count=0,
                    )

                # add to the count
                side = self.get_side(template)
                assert side is not None
                if side is Side.Left:
                    raw_sites[key].left_to_count[key.position] += 1
                elif side is Side.Right:
                    raw_sites[key].right_to_count[key.position] += 1
                else:
                    raw_sites[key].both_to_count[key.position] += 1
                raw_sites[key].count += 1

            # Aggregate sites
            all_sites: List[FindSitesMetric] = aggregate_sites(
                sites=list(raw_sites.values()), slop=inter_site_slop
            )
            # Sort the sites by reference_name
            reference_name_to_int: Dict[str, int] = {
                name: i for i, name in enumerate(reader.header.references)
            }
            all_sites = sorted(
                all_sites,
                key=lambda m: (reference_name_to_int[m.reference_name], m.position),
            )

            return all_sites

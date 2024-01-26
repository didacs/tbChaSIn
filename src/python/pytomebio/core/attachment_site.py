from attr import frozen
from samwell.dnautils import reverse_complement
from typing import ClassVar
from typing import Set


@frozen
class AttachmentSite:
    """Models an attachment site (i.e. attP and attB)

    Attributes:
        name: the name of the attachment site
        left: the sequence on the left side, in sequencing order, not including overhang
        overhang: the overhang sequence, joining the left and right
        right: the sequence on the right side, in sequencing order, not including overhang
    """

    name: str
    left: str
    overhang: str
    right: str  # typically, reverse complemented

    DNA_BASES: ClassVar[Set[str]] = {"A", "T", "G", "C"}

    @classmethod
    def parse_attachment_site(cls, site: str) -> "AttachmentSite":
        """Parse an attachment site from a string.

        The string must be a colon-delimited quadruple of
        `<name>:<left_seq>:<overhang>:<right_seq>` where `<left_seq><overhang><right_seq>` is the
        full attachment site sequence.

        Args:
            site: the string to parse
        """
        # first check that the site has >= 3 colons
        assert (
            site.count(":") == 3
        ), f"Attachment site should have four components separated by colons: name,\
                        left sequence, overhang, and right sequence.Found\
                        {site.count(':') + 1} components in: {site}"
        name, left, overhang, right = site.split(":", maxsplit=3)

        # and then check whether the string has only DNA characters remaining
        left = left.upper()
        right = right.upper()
        overhang = overhang.upper()
        assert all(
            base in cls.DNA_BASES for base in (left + overhang + right)
        ), "Attachment site has illegal bases"

        return AttachmentSite(
            name=name,
            left=left,
            overhang=overhang,
            right=reverse_complement(right),
        )


@frozen
class AttachmentSiteMatch:
    """Stores a read alignment to an attachment site.

    Attributes:
        site: the attachment site
        is_left: true if aligned to the left side of the attachment site, false for the right
        score: the alignment score
        alignment_length: the length of the alignment for the query (i.e. leading number of bases
                          to trim) including any leading query offset
    """

    site: AttachmentSite
    is_left: bool
    score: float
    alignment_length: int
    read_start_offset: int = 0

    def to_sam_tag(self) -> str:
        """Returns the SAM tag-value for this match"""
        side_label = "left" if self.is_left else "right"
        return (
            f"{self.site.name};{side_label};{self.read_start_offset};{self.alignment_length};"
            f"{self.score}"
        )

    @property
    def query_alignment_length(self) -> int:
        return self.alignment_length + self.read_start_offset

    @classmethod
    def from_sam_tag(
        cls,
        value: str,
        left: str = "",
        overhang: str = "",
        right: str = "",
    ) -> "AttachmentSiteMatch":
        """Builds an `AttachmentSiteMatch` from the given SAM tag value.

        Args:
            value: the SAM tag value to parse
            left: the sequence of the left side of the attachment site
            overhang: the overhang sequence, joining the left and right
            right: the sequence of the right side of the attachment site
        """
        try:
            name, side_label, read_start_offset, alignment_length, score = value.split(
                ";", maxsplit=4
            )
        except ValueError as ex:
            raise ValueError(f"Could not parse attachment site match SAM tag: {value}") from ex
        return AttachmentSiteMatch(
            site=AttachmentSite(name=name, left=left, overhang=overhang, right=right),
            is_left=side_label == "left",
            score=float(score),
            alignment_length=int(alignment_length),
            read_start_offset=int(read_start_offset),
        )

    """The SAM tag for the attachment site match for the current read"""
    READ_MATCH_TAG: ClassVar[str] = "rm"
    """The SAM tag for the attachment site match for the mate's read"""
    MATE_MATCH_TAG: ClassVar[str] = "mm"

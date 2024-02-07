from Bio import Align

BWA_MEM_OPEN_GAP_SCORE: int = -7
BWA_MEM_EXTEND_GAP_SCORE: int = -1


def get_global_aligner() -> Align.PairwiseAligner:
    """Aligner for aligning the full query to the full target using BWA mem inspired scoring
    parameters"""
    aligner: Align.PairwiseAligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -4
    aligner.open_gap_score = BWA_MEM_OPEN_GAP_SCORE
    aligner.extend_gap_score = BWA_MEM_EXTEND_GAP_SCORE
    return aligner


def get_glocal_aligner() -> Align.PairwiseAligner:
    """Aligner for aligning a subsequence of the query to the full target using BWA mem inspired
    scoring parameters.  Leading or trailing unaligned query bases will be represented as
    insertions.
    """
    aligner = get_global_aligner()
    aligner.mode = "local"
    aligner.query_end_open_gap_score = 0
    aligner.query_end_extend_gap_score = 0
    return aligner


def get_query_prefix_aligner() -> Align.PairwiseAligner:
    """Aligner for aligning a prefix of the query to the full target using BWA mem inspired scoring
    parameters.  Trailing unaligned query bases will be represented as insertions."""
    aligner = get_global_aligner()
    aligner.mode = "global"
    aligner.query_left_open_gap_score = BWA_MEM_OPEN_GAP_SCORE
    aligner.query_left_extend_gap_score = BWA_MEM_EXTEND_GAP_SCORE
    aligner.query_right_open_gap_score = 0
    aligner.query_right_extend_gap_score = 0
    return aligner

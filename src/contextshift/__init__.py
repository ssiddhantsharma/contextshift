"""contextshift: partition-aware functional divergence for protein families.

Give it a family of proteins and a labelling of that family from outside the
tree - genomic context, system subtype, host range - and it reports the sites
that discriminate the labels, with conservation and structural context attached.
"""

from .join import JoinThresholds, classify, summary
from .partition import Partition, Provenance, adjusted_rand_index
from .schema import ALL_SCHEMAS, SchemaError, TableSchema
from .stats import add_qvalues, benjamini_hochberg

__version__ = "0.1.0"

__all__ = [
    "ALL_SCHEMAS",
    "JoinThresholds",
    "Partition",
    "Provenance",
    "SchemaError",
    "TableSchema",
    "add_qvalues",
    "adjusted_rand_index",
    "benjamini_hochberg",
    "classify",
    "summary",
    "__version__",
]

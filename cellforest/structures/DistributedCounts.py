from dataforest.structures.DistributedContainer import DistributedContainer
from cellforest.structures.counts.Counts import Counts


class DistributedCounts(DistributedContainer):
    ELEM_CLASSES = DistributedContainer.ELEM_CLASSES + [Counts]

from dataforest.structures.DistributedContainer import DistributedContainer

from cellforest.structures.counts.Counts import Counts


class CellDistributedContainer(DistributedContainer):
    ELEM_CLASSES = DistributedContainer.ELEM_CLASSES + [Counts]

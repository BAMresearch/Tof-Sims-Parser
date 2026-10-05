from bam_masterdata.logger import logger
from bam_masterdata.metadata.entities import CollectionType


class TestOpenbisParserExample:
    def test_parse(self, parser):
        collection = CollectionType()
        parser.parse([], collection, logger)

        assert len(collection.attached_objects) == 0
        objects = list(collection.attached_objects.values())
        objects
        assert len(collection.relationships) == 0

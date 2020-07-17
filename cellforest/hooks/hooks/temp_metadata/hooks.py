import os

from dataforest.hooks import hook

from cellforest.templates.WriterMethodsSC import WriterMethodsSC


@hook(attrs=["temp_meta"])
def hook_store_temp_meta(dp):
    """
    Stores a temporary metadata file for the current process run, which will be
    removed by `hook_clean_temp_metadata`
    """
    if dp.temp_meta:
        dp._metadata_filepath = dp.forest[dp.process_name].path / dp.forest.schema.__class__.TEMP_METADATA_FILENAME
        WriterMethodsSC.tsv(
            dp._metadata_filepath, dp.forest[dp.process_name].forest.meta, header=True,
        )


@hook(attrs=["temp_meta"])
def hook_clean_temp_meta(dp):
    """Removes temporary metadata file for the current process run"""
    if dp.temp_meta:
        os.remove(dp._metadata_filepath)

import os

from dataforest.hooks import hook, dataprocess

from cellforest.templates.WriterMethodsSC import WriterMethodsSC


@hook(attrs=["temp_meta"])
def hook_store_temp_meta(dp: dataprocess):
    """
    Stores a temporary metadata file for the current process run, which will be
    removed by `hook_clean_temp_metadata`
    """
    if dp.temp_meta:
        dp._metadata_filepath = dp.branch[dp.name].path / dp.branch.schema.__class__.TEMP_METADATA_FILENAME
        process_name = dp.branch.current_process
        precursors = dp.branch.spec.get_precursors_lookup()[process_name]
        precursor = precursors[-1] if precursors else None
        temp_meta = dp.branch[dp.name].branch._get_meta(precursor)
        temp_meta = dp.branch._apply_data_ops(dp.name, temp_meta)
        WriterMethodsSC.tsv(dp._metadata_filepath, temp_meta, header=True)


@hook(attrs=["temp_meta"])
def hook_clean_temp_meta(dp: dataprocess):
    """Removes temporary metadata file for the current process run"""
    if dp.temp_meta:
        os.remove(dp._metadata_filepath)


@hook
def hook_clear_metadata_cache(dp: dataprocess):
    """
    Clear metadata attribute after process run so that it can be recalculated
    """
    dp.branch._meta = None

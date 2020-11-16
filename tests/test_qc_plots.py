import pytest
import os
from PIL import Image

from tests.fixtures import *
import cellforest as cf


def test_root_plots(load_test_config, remove_plots, build_root):
    plot_map = build_root["root"].plot_map
    assert "_UMIS_PER_CELL_HIST_" in plot_map
    assert len(plot_map.keys()) == 1


def test_pyplot_plots(load_test_config, remove_plots, test_norm):
    plot_map = test_norm["normalize"].plot_map
    plot_lookup = test_norm["normalize"]._plot_lookup
    plot_name = "_GENES_PER_CELL_HIST_"

    assert plot_name in plot_lookup
    assert "genes_per_cell_hist-plot_size:800+800-stratify:sample.png" in plot_lookup[plot_name].values()  # valid
    assert "genes_per_cell_hist-plot_size:1600+800-stratify:none.png" in plot_lookup[plot_name].values()  # valid


def test_ggplot2_plots(load_test_config, remove_plots, test_cluster):
    plot_map = test_cluster["cluster"].plot_map
    plot_lookup = test_cluster["cluster"]._plot_lookup
    plot_name = "_PERC_HSP_PER_CELL_VLN"

    assert plot_name in plot_lookup
    assert (
        "perc_hsp_per_cell_vln-plot_size:1600+800-stratify:cluster_id.png" in plot_lookup[plot_name].values()
    )  # valid

    assert os.listdir(plot_map[plot_name]['{"plot_size": "wide", "stratify": "cluster"}'].parent) == [
        "perc_hsp_per_cell_vln-plot_size:1600+800-stratify:cluster_id.png"
    ]


def test_infer_plot_names(load_test_config, remove_plots, build_root):
    plot_map = build_root["root"].plot_map
    plot_lookup = build_root["root"]._plot_lookup
    plot_name = "_UMIS_PER_CELL_HIST_"

    assert plot_name in plot_map
    assert plot_name in plot_lookup
    assert '{"plot_size": "default", "stratify": "default"}' in plot_map[plot_name]
    assert '{"plot_size": "default", "stratify": "default"}' in plot_lookup[plot_name]

    # check that the inferred name is correct
    assert (
        plot_lookup[plot_name]['{"plot_size": "default", "stratify": "default"}']
        == "umis_per_cell_hist-plot_size:800+800-stratify:none.png"
    )

    # check that only one correct plot is in the _plots directory
    plot_dir = list(plot_map[plot_name].values())[0].parent
    assert os.listdir(plot_dir) == ["umis_per_cell_hist-plot_size:800+800-stratify:none.png"]


def test_invalid_plot(load_test_config, remove_plots, test_norm):
    plot_map = test_norm["normalize"].plot_map
    plot_lookup = test_norm["normalize"]._plot_lookup
    plot_name = "_GENES_PER_CELL_HIST_"

    assert plot_name in plot_lookup
    assert (
        "genes_per_cell_hist-non_existent.png" in plot_lookup[plot_name].values()
    )  # TODO-QC: this file should ideally be removed from plot_map
    # but, this file should not be in plot_map
    assert sorted(
        os.listdir(plot_map[plot_name]['{"plot_size": "large", "stratify": "nonexistent_column"}'].parent)
    ) == sorted(
        [
            "genes_per_cell_hist-plot_size:800+800-stratify:sample.png",
            "genes_per_cell_hist-plot_size:1600+800-stratify:none.png",
        ]
    )


def test_additional_kwargs(load_test_config, remove_plots, test_reduce):
    plot_map = test_reduce["reduce"].plot_map
    plot_lookup = test_reduce["reduce"]._plot_lookup
    plot_name = "_PCA_EMBEDDINGS_SCAT_"

    assert '{"alpha": 0.4, "npcs": 7, "plot_size": "large", "stratify": "nFeature_RNA"}' in plot_map[plot_name]
    assert '{"alpha": 0.4, "npcs": 7, "plot_size": "large", "stratify": "entity_id"}' in plot_map[plot_name]
    assert '{"alpha": 0.4, "npcs": 7, "plot_size": "large", "stratify": "nFeature_RNA"}' in plot_lookup[plot_name]
    assert '{"alpha": 0.4, "npcs": 7, "plot_size": "large", "stratify": "entity_id"}' in plot_lookup[plot_name]

    assert (
        "pca_embeddings_scat-alpha:0.4-npcs:7-plot_size:1600+1600-stratify:nfeature_rna.png"
        in plot_lookup[plot_name].values()
    )
    assert (
        "pca_embeddings_scat-alpha:0.4-npcs:7-plot_size:1600+1600-stratify:sample.png"
        in plot_lookup[plot_name].values()
    )

    plot_dir = list(plot_map[plot_name].values())[0].parent
    assert sorted(os.listdir(plot_dir)) == sorted(
        [
            "pca_embeddings_scat-alpha:0.4-npcs:7-plot_size:1600+1600-stratify:nfeature_rna.png",
            "pca_embeddings_scat-alpha:0.4-npcs:7-plot_size:1600+1600-stratify:sample.png",
        ]
    )


def test_plot_size(load_test_config, remove_plots, build_root):
    plot_map = build_root["root"].plot_map
    assert "_UMIS_PER_CELL_HIST_" in plot_map

    path_to_plot_file = list(plot_map["_UMIS_PER_CELL_HIST_"].values())[0]
    im = Image.open(path_to_plot_file)
    width, height = im.size

    assert width == 800
    assert height == 800

from __future__ import annotations

from pathlib import Path

from matplotlib.patches import Patch

from phytreeviz import TreeViz, load_example_tree_file


def test_treeviz1(tmp_path: Path):
    """Test TreeViz"""
    outfile = tmp_path / "result.png"

    tree_file = load_example_tree_file("small_example.nwk")

    tv = TreeViz(tree_file)
    tv.show_branch_length(color="red")
    tv.show_confidence(color="blue")
    tv.show_scale_bar()

    tv.savefig(outfile, dpi=300)
    assert outfile.exists()


def test_treeviz2(tmp_path: Path):
    """Test TreeViz"""
    outfile = tmp_path / "result.png"

    tree_file = load_example_tree_file("small_example.nwk")

    tv = TreeViz(tree_file, height=0.7)
    tv.show_scale_axis()

    tv.set_node_label_props("Homo_sapiens", color="grey")
    tv.set_node_label_props("Pongo_abelii", color="green", style="italic")

    tv.set_node_line_props(
        ["Hylobates_moloch", "Nomascus_leucogenys"], color="orange", lw=2
    )
    tv.set_node_line_props(
        ["Homo_sapiens", "Pan_troglodytes", "Pan_paniscus"],
        color="magenta",
        ls="dotted",
    )

    tv.savefig(outfile, dpi=300)
    assert outfile.exists()


def test_treeviz3(tmp_path: Path):
    """Test TreeViz"""
    outfile = tmp_path / "result.png"

    tree_file = load_example_tree_file("small_example.nwk")

    tv = TreeViz(tree_file, align_leaf_label=True)
    tv.show_scale_axis()

    group1 = ["Hylobates_moloch", "Nomascus_leucogenys"]
    group2 = ["Homo_sapiens", "Pan_paniscus"]

    tv.highlight(group1, "orange")
    tv.highlight(group2, "lime")

    tv.annotate(group1, "group1")
    tv.annotate(group2, "group2")

    tv.marker(group1, marker="s", color="blue")
    tv.marker(group2, marker="D", color="purple", descendent=True)
    tv.marker("Pongo_abelii", color="red")

    tv.savefig(outfile, dpi=300)
    assert outfile.exists()


def test_treeviz4(tmp_path: Path):
    """Test TreeViz"""
    outfile = tmp_path / "result.png"

    tree_file = load_example_tree_file("medium_example.nwk")

    tv = TreeViz(tree_file, height=0.3, align_leaf_label=True, leaf_label_size=10)
    tv.show_scale_bar()

    group1 = ["Hylobates_moloch", "Nomascus_leucogenys"]
    group2 = ["Homo_sapiens", "Pongo_abelii"]
    group3 = ["Piliocolobus_tephrosceles", "Rhinopithecus_bieti"]
    group4 = ["Chlorocebus_sabaeus", "Papio_anubis"]

    tv.highlight(group1, "orange", area="full")
    tv.highlight(group2, "skyblue", area="full")
    tv.highlight(group3, "lime", area="full")
    tv.highlight(group4, "pink", area="full")

    tv.link(group3, group4, connectionstyle="arc3,rad=0.2")

    fig = tv.plotfig()

    _ = fig.legend(
        handles=[
            Patch(label="group1", color="orange"),
            Patch(label="group2", color="skyblue"),
            Patch(label="group3", color="lime"),
            Patch(label="group4", color="pink"),
        ],
        frameon=False,
        bbox_to_anchor=(0.3, 0.3),
        loc="center",
        ncols=2,
    )

    fig.savefig(str(outfile), dpi=300)
    assert outfile.exists()

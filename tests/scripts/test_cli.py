from __future__ import annotations

import subprocess as sp
from pathlib import Path

from phytreeviz import load_example_tree_file


def test_phytreeviz_cli1(tmp_path: Path):
    """Test phyTreeViz CLI"""
    outfile = tmp_path / "result.png"

    cmd = f"phytreeviz -i '((A,B),((C,D),(E,(F,G))));' -o {outfile}"
    sp.run(cmd, shell=True)

    assert outfile.exists()


def test_phytreeviz_cli2(tmp_path: Path):
    """Test phyTreeViz CLI"""
    outfile = tmp_path / "result.png"

    tree_file = load_example_tree_file("small_example.nwk")
    cmd = f"phytreeviz -i {tree_file} -o {outfile} "
    cmd += "--show_branch_length --show_confidence "
    sp.run(cmd, shell=True)

    assert outfile.exists()


def test_phytreeviz_cli3(tmp_path: Path):
    """Test phyTreeViz CLI"""
    outfile = tmp_path / "result.png"

    tree_file = load_example_tree_file("medium_example.nwk")
    cmd = f"phytreeviz -i {tree_file} -o {outfile} "
    cmd += "--fig_height 0.3 --align_leaf_label "
    sp.run(cmd, shell=True)

    assert outfile.exists()

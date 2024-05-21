from pathlib import Path
from typing import List
from .molecule import Molecule
import operator
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Conversion of the second argument of issubdtype")


def write_reads_or_molecules(path, mols_or_reads, require_umis=True, sort=True):
    with open(path, "w") as f:
        if require_umis:
            if sort:
                mols_or_reads = sorted(
                    mols_or_reads,
                    key=lambda mol_or_read: (
                        mol_or_read.umi,
                        mol_or_read.cell_id,
                        mol_or_read.clone_id,
                    ),
                )
            print("#cell_id", "umi", "clone_id", sep="\t", file=f)
            for mol_or_read in mols_or_reads:
                print(
                    mol_or_read.cell_id,
                    mol_or_read.umi,
                    mol_or_read.clone_id,
                    sep="\t",
                    file=f,
                )
        else:
            if sort:
                mols_or_reads = sorted(
                    mols_or_reads,
                    key=lambda mol_or_read: (mol_or_read.clone_id, mol_or_read.cell_id),
                )
            print("#cell_id", "clone_id", sep="\t", file=f)
            for mol_or_read in mols_or_reads:
                print(mol_or_read.cell_id, mol_or_read.clone_id, sep="\t", file=f)


def write_adata(molecules: List[Molecule], output):
    """
    Creates Annotated Data matrix object (.h5ad) with UMI count matrix.
    """
    from scipy.sparse import csr_matrix
    import anndata as ad
    import pandas as pd
    import numpy as np

    # Creation of indexes in sparse matrix
    cell_ids = list(set([molecule.cell_id for molecule in molecules]))
    clone_ids = list(set([molecule.clone_id for molecule in molecules]))
    cell_idx = dict(zip(cell_ids, range(len(cell_ids))))
    clone_idx = dict(zip(clone_ids, range(len(clone_ids))))

    # Filling the sparse matrix
    matrix = csr_matrix(np.zeros((len(cell_ids), len(clone_ids)), dtype=int))
    for molecule in molecules:
        matrix[cell_idx[molecule.cell_id], clone_idx[molecule.clone_id]] += 1

    ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=cell_ids),
        var=pd.DataFrame(index=clone_ids),
        dtype=int,
    ).write_h5ad(output)
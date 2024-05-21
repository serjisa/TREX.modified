"""
Run on 10X data
"""
import sys
import logging
from pathlib import Path
from collections import Counter
from typing import List, Dict, Iterable

import pandas as pd

from . import (
    setup_logging,
    CommandLineError,
    add_file_logging,
    make_output_dir,
    add_common_arguments,
)
from .. import __version__
from ..cellranger import make_cellranger, CellRangerError
from ..writers import (
    write_reads_or_molecules,
    write_adata,
)
from ..clustering import cluster_sequences
from ..molecule import Molecule, compute_molecules
from ..error import TrexError
from ..dataset import DatasetReader


__author__ = "leonie.von.berlin@ki.se"
# modifications: sergey.isaev@meduniwien.ac.at

logger = logging.getLogger(__name__)


def main(args):
    setup_logging(debug=args.debug)

    output_dir = args.output
    try:
        make_output_dir(output_dir, args.delete)
    except FileExistsError:
        raise CommandLineError(
            f"Output directory '{output_dir}' already exists "
            "(use --delete to force deleting an existing output directory)"
        )

    add_file_logging(output_dir / "log.txt")
    logger.info(f"Trex {__version__}")
    logger.info("Command line arguments: %s", " ".join(sys.argv[1:]))

    allowed_cell_ids = None
    if args.filter_cellids:
        allowed_cell_ids = read_allowed_cellids(args.filter_cellids)
    transcriptome_inputs = args.path
    if args.samples:
        sample_names = args.samples.split(",")
    elif len(transcriptome_inputs) == 1:
        sample_names = [None]  # Do not modify suffixes
    else:
        sample_names = [path.name for path in transcriptome_inputs]
        logger.info("Using these sample names: %s", ", ".join(sample_names))
    if len(sample_names) != len(transcriptome_inputs):
        raise CommandLineError(
            "The number of sample names (--samples) must match the number of "
            "provided transcriptome datasets"
        )
    if args.amplicon:
        amplicon_inputs = args.amplicon
        if len(transcriptome_inputs) != len(amplicon_inputs):
            raise CommandLineError(
                "As many amplicon as transcriptome datasets must be provided"
            )
    else:
        amplicon_inputs = []

    try:
        run_trex(
            output_dir,
            genome_name=args.genome_name,
            allowed_cell_ids=allowed_cell_ids,
            chromosome=args.chromosome,
            start=args.start - 1 if args.start is not None else None,
            end=args.end,
            transcriptome_inputs=transcriptome_inputs,
            amplicon_inputs=amplicon_inputs,
            sample_names=sample_names,
            prefix=args.prefix,
            max_hamming=args.max_hamming,
            max_dashes=args.max_dashes,
        )
    except (CellRangerError, TrexError) as e:
        raise CommandLineError(e)


def add_arguments(parser):
    groups = add_common_arguments(parser)


def run_trex(
    output_dir: Path,
    genome_name: str,
    allowed_cell_ids: List[str],
    chromosome: str,
    start: int,
    end: int,
    transcriptome_inputs: List[Path],
    amplicon_inputs: List[Path],
    sample_names: List[str],
    prefix: bool,
    max_hamming: int,
    max_dashes: int,
):
    if len(sample_names) != len(set(sample_names)):
        raise TrexError("The sample names need to be unique")

    dataset_reader = DatasetReader(
        output_dir, genome_name, chromosome, start, end, prefix
    )
    reads = dataset_reader.read_all(
        transcriptome_inputs, amplicon_inputs, sample_names, allowed_cell_ids
    )
    if not reads:
        raise TrexError("No reads left after --filter-cellids filtering")

    clone_ids = [
        r.clone_id for r in reads if "-" not in r.clone_id and "0" not in r.clone_id
    ]
    logger.info(
        f"Read {len(reads)} reads containing (parts of) the cloneID "
        f"({len(clone_ids)} full cloneIDs, {len(set(clone_ids))} unique)"
    )

    write_reads_or_molecules(output_dir / "reads.txt", reads)

    molecules = compute_molecules(reads)
    clone_ids = [
        m.clone_id for m in molecules if "-" not in m.clone_id and "0" not in m.clone_id
    ]
    logger.info(
        f"Detected {len(molecules)} molecules ({len(clone_ids)} full cloneIDs, "
        f"{len(set(clone_ids))} unique)"
    )

    write_reads_or_molecules(output_dir / "molecules.txt", molecules, sort=False)
    
    corrected_molecules = correct_clone_ids_per_cell(molecules, max_hamming, max_dashes)
    corrected_molecules = correct_gaps(corrected_molecules)
    
    clone_ids = [
        m.clone_id
        for m in corrected_molecules
        if "-" not in m.clone_id and "0" not in m.clone_id
    ]
    logger.info(
        f"After cloneID correction, {len(set(clone_ids))} unique cloneIDs remain"
    )
    
    logger.info("Writing UMI matrix")
    write_adata(corrected_molecules, output_dir / "umi_matrix.h5ad")
    
    
def correct_gaps(molecules: List[Molecule]) -> List[Molecule]:
    """
    Correction of CloneID molecules with gaps. The function corrects CloneID
    sequence with gaps if there's only one sequence without gaps that matches with
    target sequence.
    """
    clone_ids = list(set([m.clone_id for m in molecules]))
    clone_ids_gaps = [clone_id for clone_id in clone_ids if "-" in clone_id]
    clone_ids_wo_gaps = pd.Series([clone_id for clone_id in clone_ids if "-" not in clone_id])
    
    corrected_clone_ids = {}
    blacklist = []
    for clone_id_gap in clone_ids_gaps:
        if clone_id_gap[0] == "-":
            mathed_clone_ids = list(
                clone_ids_wo_gaps[clone_ids_wo_gaps.str.endswith(clone_id_gap.lstrip("-"))]
            )
            if len(mathed_clone_ids) == 1:
                corrected_clone_ids[clone_id_gap] = mathed_clone_ids[0]
            elif len(mathed_clone_ids) == 0:
                corrected_clone_ids[clone_id_gap] = clone_id_gap
            else:
                blacklist.append(clone_id_gap)
        elif clone_id_gap[-1] == "-":
            mathed_clone_ids = list(
                clone_ids_wo_gaps[clone_ids_wo_gaps.str.startswith(clone_id_gap.rstrip("-"))]
            )
            if len(mathed_clone_ids) == 1:
                corrected_clone_ids[clone_id_gap] = mathed_clone_ids[0]
            elif len(mathed_clone_ids) == 0:
                corrected_clone_ids[clone_id_gap] = clone_id_gap
            else:
                blacklist.append(clone_id_gap)
        else:
            corrected_clone_ids[clone_id_gap] = clone_id_gap
            
    corrected_molecules = []
    for molecule in molecules:
        if molecule.clone_id not in blacklist:
            if molecule.clone_id in corrected_clone_ids:
                corrected_molecules.append(
                    molecule._replace(clone_id=corrected_clone_ids[molecule.clone_id])
                )
            else:
                corrected_molecules.append(molecule)
    return corrected_molecules


def read_allowed_cellids(path):
    """
    Read a user-provided list of allowed cell IDs from a CSV

    Example:

    "X","z"
    1,"ACGTACGTACGTACGT_10x99"

    or:

    "X","z"
    1,"ACGTACGTACGTACGT"
    """
    allowed_ids = []
    filtered_df = pd.read_csv(Path(path), sep=",", index_col=0)
    for cell_id in filtered_df.iloc[:, 0]:
        if cell_id.endswith("-1"):
            raise TrexError(
                "Cell ids in the list of allowed cell IDs must not end in '-1'"
            )
        allowed_ids.append(cell_id)
    logger.info(f"Restricting analysis to {len(allowed_ids)} allowed cells")
    return set(allowed_ids)


def correct_clone_ids(
    molecules: List[Molecule], max_hamming: int
) -> List[Molecule]:
    """
    Attempt to correct sequencing errors in the cloneID sequences of all molecules.
    """
    import numpy as np
    
    # Obtain all cloneIDs (including those with '-' and '0')
    clone_ids = [m.clone_id for m in molecules]

    def hamming(s1, s2):
        # Manual hamming distance function: fast and works dashes and 0
        h = 0
        for c1, c2 in zip(s1, s2):
            if (c1 != c2) & (c1 != "-") & (c2 != "-"):
                h += 1
        return h

    # Similarity function
    def is_similar(s1, s2):
        n_dash = sum([ ((c1 == "-") or (c2 == "-")) for (c1, c2) in zip(s1, s2) ])
        # Current hamming is a Hamming distance with adjustment for number of "-"
        current_max_hamming = max_hamming - n_dash * max_hamming / 30
        return hamming(s1, s2) <= current_max_hamming

    clusters = cluster_sequences(list(set(clone_ids)), is_similar=is_similar)
    
    # Getting consensus for cluster of sequences
    def consensus(s, weight_nonATGC=0.1):
        # Consensus instead of the most representative
        c = ""
        for i in range(len(s[0])):
            letters = pd.Series([seq[i] for seq in s]).value_counts()
            if "-" in letters.index:
                letters["-"] = letters["-"] * weight_nonATGC
            if "0" in letters.index:
                letters["0"] = letters["0"] * weight_nonATGC
            c += letters.sort_values().index[-1]
        return c

    # Map non-singleton cloneIDs to a cluster representative
    clone_id_map = dict()
    for cluster in clusters:
        if len(cluster) > 1:
            # Pick consensus of cloneIDs as representative. For this reason get
            # all of the clone IDs in cluster with duplications
            whole_cluster_mask = np.isin(clone_ids, cluster)
            whole_cluster = [
                clone_ids[i] for i in range(len(clone_ids)) if whole_cluster_mask[i]
            ]
            representative = consensus(whole_cluster)
            for clone_id in cluster:
                clone_id_map[clone_id] = representative

    # Create a new list of molecules in which the cloneIDs have been replaced
    # by consensus within cluster
    new_molecules = []
    for molecule in molecules:
        clone_id = clone_id_map.get(molecule.clone_id, molecule.clone_id)
        new_molecules.append(molecule._replace(clone_id=clone_id))
        
    return new_molecules


def correct_clone_ids_per_cell(
    molecules: List[Molecule], max_hamming: int, max_dashes: int = 10 
) -> List[Molecule]:
    """
    Performs cloneID correction in each cell independently. It's better than
    whole pseudo-transcriptome correction because of two reasons:
    
    (1) with increasing of number of samples there will be artificial bridges
    between cloneID clusters,
    (2) usually set of cloneIDs within cell is good enough for such a correction.
    """
    from collections import defaultdict

    # Creation of a dictionary with molecules per cell
    cells = defaultdict(list)
    for molecule in molecules:
        if molecule.clone_id.count("-") <= max_dashes:
            cells[molecule.cell_id].append(molecule)

    # Correction of molecules within each cell
    molecules_corrected = []
    for cell_id in cells:
        molecules_corrected += correct_clone_ids(
            cells[cell_id],
            max_hamming=max_hamming,
        )
    
    return molecules_corrected
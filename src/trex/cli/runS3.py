"""
Run on SmartSeq3 data
"""
import sys
import operator
import shutil
import warnings
import logging
from pathlib import Path
from collections import Counter
from typing import List, Dict, Iterable

from tinyalign import hamming_distance
import numpy as np
import pandas as pd
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Conversion of the second argument of issubdtype')
    import loompy

from .. import __version__
from ..utils import NiceFormatter
from ..clustering import cluster_sequences
from ..clone import CloneGraph
from ..cell import Cell, compute_cells
from ..error import BraintraceError
from ..dataset import DatasetReader
from ..bam import Read


__author__ = 'leonie.von.berlin@ki.se'

logger = logging.getLogger(__name__)


def main(args):
    setup_logging(debug=args.debug)

    output_dir = args.output
    try:
        make_output_dir(output_dir, args.delete)
    except FileExistsError:
        logger.error(f'Output directory "{output_dir}" already exists '
                     '(use --delete to force deleting an existing output directory)')
        sys.exit(1)

    add_file_logging(output_dir / 'log.txt')
    logger.info(f'Braintrace {__version__}')
    logger.info('Command line arguments: %s', ' '.join(sys.argv[1:]))

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
        logger.error("The number of sample names (--samples) must match the number of "
            "provided transcriptome datasets")
        sys.exit(1)
    if args.amplicon:
        amplicon_inputs = args.amplicon
        if len(transcriptome_inputs) != len(amplicon_inputs):
            logger.error("As many amplicon as transcriptome datasets must be provided")
            sys.exit(1)
    else:
        amplicon_inputs = []

    highlight_cell_ids = []
    if args.highlight:
        with open(args.highlight) as f:
            highlight_cell_ids = [line.strip() for line in f]
    
    try:
        run_braintrace(
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
            min_length=args.min_length,
            jaccard_threshold=args.jaccard_threshold,
            keep_single_reads=args.keep_single_reads,
            should_write_read_matrix=args.read_matrix,
            should_plot=args.plot,
            highlight_cell_ids=highlight_cell_ids,
        )
    except (BraintraceError) as e:
        logger.error("%s", e)
        sys.exit(1)


def add_arguments(parser):
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--debug', default=False, action='store_true',
        help='Print some extra debugging messages')
    parser.add_argument('--genome-name', metavar='NAME',
        help='Name of the genome as indicated in Cell Ranger count run with the flag --genome. '
             'Default: None',
        default=None)
    parser.add_argument('--chromosome', '--chr',
        help='Name of chromosome on which clone ID is located.'
             ' Default: Last chromosome in BAM file',
        default=None)
    parser.add_argument('--output', '-o', '--name', '-n', metavar='DIRECTORY', type=Path,
        help='name of the run and directory created by program. Default: %(default)s',
        default=Path('braintrace_run'))
    parser.add_argument('--delete', action='store_true',
        help='Delete output directory if it exists')
    parser.add_argument('--start', '-s',
        help='Position of first clone ID nucleotide. Default: Auto-detected',
        type=int, default=None)
    parser.add_argument('--end', '-e',
        help='Position of last clone ID nucleotide. Default: Auto-detected',
        type=int, default=None)
    parser.add_argument('--min-length', '-m',
        help='Minimum number of nucleotides a clone ID must have. Default: %(default)s',
        type=int, default=20)
    parser.add_argument('--max-hamming',
        help='Maximum hamming distance allowed for two clone IDs to be called similar. '
            'Default: %(default)s',
        type=int, default=5)
    parser.add_argument('--jaccard-threshold', type=float, default=0, metavar='VALUE',
        help='If the Jaccard index between clone IDs of two cells is higher than VALUE, the cells '
            'are considered similar. Default: %(default)s')
    parser.add_argument('--amplicon', '-a', nargs='+', metavar='DIRECTORY',
        help='Path to united bam file for all cells or path to a folder with one bam file per cell,'
            'containing sequencing of the clone ID amplicon library. Provide these in '
            'same order as transcriptome datasets', default=None)
    parser.add_argument('--filter-cellids', '-f', metavar='CSV', type=Path,
        help='CSV file containing cell IDs to keep in the analysis.'
             ' This flag enables to remove cells e.g. doublets',
        default=None)
    parser.add_argument('--keep-single-reads', action='store_true', default=False,
        help='Keep clone IDs supported by only a single read. Default: Discard them')
    parser.add_argument('--highlight',
        help='Highlight cell IDs listed in FILE in the clone graph')
    parser.add_argument('--samples',
        help='Sample names separated by comma, in the same order as bam files/bam file directories',
        default=None)
    parser.add_argument("--prefix", default=False, action="store_true",
        help="Add sample name as prefix to cell IDs (instead of as suffix)")
    parser.add_argument('--read-matrix', default=False, action='store_true',
        help='Creates a read count matrix with cells as columns and clone IDs as rows')
    parser.add_argument('--plot', dest='plot', default=False, action='store_true',
        help='Plot the clone graph')
    parser.add_argument('path', type=Path, nargs='+', metavar='DIRECTORY',
        help='Path to a united bam file for all cells or path to a folder with one bam file per cell.')


def setup_logging(debug: bool) -> None:
    """
    Set up logging. If debug is True, then DEBUG level messages are printed.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(NiceFormatter())

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)


def add_file_logging(path: Path) -> None:
    file_handler = logging.FileHandler(path)
    root = logging.getLogger()
    root.addHandler(file_handler)


def make_output_dir(path, delete_if_exists):
    try:
        path.mkdir()
    except FileExistsError:
        if delete_if_exists:
            logger.debug(f'Re-creating folder "{path}"')
            shutil.rmtree(path)
            path.mkdir()
        else:
            raise


def run_braintrace(
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
    min_length: int,
    jaccard_threshold: float,
    keep_single_reads: bool,
    should_write_read_matrix: bool,
    should_plot: bool,
    highlight_cell_ids: List[str],
):
    if len(sample_names) != len(set(sample_names)):
        raise BraintraceError("The sample names need to be unique")

    dataset_reader = DatasetReader(output_dir, genome_name, chromosome, start, end, prefix)
    reads = dataset_reader.read_all(
        transcriptome_inputs, amplicon_inputs, sample_names, allowed_cell_ids, require_umis=False, cell_id_tag = "BC")
    if not reads:
        raise BraintraceError("No reads left after --filter-cellids filtering")

    clone_ids = [
        r.clone_id for r in reads if '-' not in r.clone_id and '0' not in r.clone_id]
    logger.info(f'Read {len(reads)} reads containing (parts of) the clone ID '
        f'({len(clone_ids)} full clone IDs, {len(set(clone_ids))} unique)')

    write_reads(output_dir / "reads.txt", reads)

    corrected_reads = correct_clone_ids(reads, max_hamming, min_length)
    clone_ids = [r.clone_id for r in corrected_reads
        if '-' not in r.clone_id and '0' not in r.clone_id]
    logger.info(f'After clone ID correction, {len(set(clone_ids))} unique clone IDs remain')

    write_reads(output_dir / 'reads_corrected.txt', corrected_reads)

    # We do not have multiple reads per molecule, so we treat each read as one molecule
    corrected_molecules = corrected_reads
    cells = compute_cells(corrected_molecules, min_length)
    logger.info(f'Detected {len(cells)} cells')
    write_cells(output_dir / 'cells.txt', cells)

    #cells = filter_cells(cells, corrected_reads, keep_single_reads)
    #logger.info(f'{len(cells)} filtered cells remain')
    #write_cells(output_dir / 'cells_filtered.txt', cells)

    if should_write_read_matrix:
        logger.info(f"Writing read matrix")
        write_read_matrix(output_dir, cells)

    clone_graph = CloneGraph(cells, jaccard_threshold=jaccard_threshold)

    with open(output_dir / 'components.txt', 'w') as components_file:
        print(clone_graph.components_txt(highlight_cell_ids), file=components_file, end='')
    if should_plot:
        logger.info('Plotting clone graph')
        clone_graph.plot(output_dir / 'graph', highlight_cell_ids)

    bridges = clone_graph.bridges()
    logger.info(f'Removing {len(bridges)} bridges from the graph')
    clone_graph.remove_edges(bridges)
    with open(output_dir / 'components_corrected.txt', 'w') as components_file:
        print(clone_graph.components_txt(highlight_cell_ids), file=components_file, end='')

    if should_plot:
        logger.info('Plotting corrected clone graph')
        clone_graph.plot(output_dir / 'graph_corrected', highlight_cell_ids)

    clones = clone_graph.clones()
    with open(output_dir / 'clones.txt', 'w') as f:
        clone_graph.write_clones(f, clones)
    with open(output_dir / 'clone_sequences.txt', 'w') as f:
        clone_graph.write_clone_sequences(f, clones)
    logger.info(f'Detected {len(clones)} clones')
    clone_sizes = Counter(len(cells) for clone_id, cells in clones)
    logger.info('Clone size histogram\n size count\n%s',
        '\n'.join(f'{k:5d} {clone_sizes[k]:5d}' for k in sorted(clone_sizes)))
    number_of_cells_in_clones = sum(k * v for k, v in clone_sizes.items())
    logger.debug('No. of cells in clones: %d', number_of_cells_in_clones)
    assert len(cells) == number_of_cells_in_clones



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
            raise BraintraceError("Cell ids in the list of allowed cell IDs must not end in '-1'")
        allowed_ids.append(cell_id)
    logger.info(f'Restricting analysis to {len(allowed_ids)} allowed cells')
    return set(allowed_ids)

#EDIT
def correct_clone_ids(
        reads: List[Read], max_hamming: int, min_overlap: int = 20) -> List[Read]:
    """
    Attempt to correct sequencing errors in the clone ID sequences of all molecules
    """
    # Obtain all clone IDs (including those with '-' and '0')
    clone_ids = [r.clone_id for r in reads]

    # Count the full-length clone IDs
    clone_id_counts = Counter(clone_ids)

    # Cluster them by Hamming distance
    def is_similar(s, t):
        # m = max_hamming
        if '-' in s or '-' in t:
            # Remove suffix and/or prefix where sequences do not overlap
            s = s.lstrip('-')
            t = t[-len(s):]
            s = s.rstrip('-')
            if len(s) < min_overlap:
                return False
            t = t[:len(s)]
            # TODO allowed Hamming distance should be reduced relative to the overlap length
            # m = max_hamming * len(s) / len(original_length_of_s)
        return hamming_distance(s, t) <= max_hamming

    clusters = cluster_sequences(list(set(clone_ids)), is_similar=is_similar, k=7)

    # Map non-singleton clone IDs to a cluster representative
    clone_id_map = dict()
    for cluster in clusters:
        if len(cluster) > 1:
            # Pick most frequent clone ID as representative
            representative = max(cluster, key=lambda bc: (clone_id_counts[bc], bc))
            for clone_id in cluster:
                clone_id_map[clone_id] = representative

    # Create a new list of molecules in which the clone IDs have been replaced
    # by their representatives
    new_reads = []
    for read in reads:
        clone_id = clone_id_map.get(read.clone_id, read.clone_id)
        new_reads.append(read._replace(clone_id=clone_id))
    return new_reads
'''
#EDIT
def filter_cells(
    cells: Iterable[Cell],
    reads: Iterable[Read],
    keep_single_reads: bool = False
) -> List[Cell]:
    """
    Filter clone IDs according to two criteria:

    - Clone ids that have only a count of one and can be found in another cell are most
      likely results of contamination and are removed,
    - If keep_single_reads is False, clone IDs that have only a count of one and are also only based
      on one read are also removed
    """
    overall_clone_id_counts: Dict[str, int] = Counter()
    for cell in cells:
        overall_clone_id_counts.update(cell.clone_id_counts)

    single_read_clone_ids = set()
    for molecule in molecules:
        if molecule.read_count == 1:
            single_read_clone_ids.add(molecule.clone_id)
    logger.info(f"Found {len(single_read_clone_ids)} single-read clone IDs")

    # filters out clone IDs with a count of one that appear in another cell
    new_cells = []
    for cell in cells:
        clone_id_counts = cell.clone_id_counts.copy()
        for clone_id, count in cell.clone_id_counts.items():
            if count > 1:
                # This clone ID occurs more than once in this cell - keep it
                continue
            if overall_clone_id_counts[clone_id] > 1:
                # This clone ID occurs also in other cells - remove it
                del clone_id_counts[clone_id]
            elif clone_id in single_read_clone_ids and not keep_single_reads:
                del clone_id_counts[clone_id]
        if clone_id_counts:
            new_cells.append(Cell(cell_id=cell.cell_id, clone_id_counts=clone_id_counts))
    return new_cells
'''

def write_reads(path, reads):
    with open(path, 'w') as f:
        print("#cell_id", "clone_id", sep="\t", file=f)
        for read in sorted(reads, key=lambda read: (read.clone_id, read.cell_id)):
            print(read.cell_id, read.clone_id, sep='\t', file=f)


def write_cells(path: Path, cells: List[Cell]) -> None:
    """Write cells to a tab-separated file"""
    with open(path, 'w') as f:
        print(
            "#cell_id", ":", "clone_id1", "count1", "clone_id2", "count2", "...", sep="\t", file=f)
        for cell in cells:
            row = [cell.cell_id, ':']
            sorted_clone_ids = sorted(
                cell.counts, key=lambda x: cell.counts[x], reverse=True)
            if not sorted_clone_ids:
                continue
            for clone_id in sorted_clone_ids:
                row.extend([clone_id, cell.counts[clone_id]])
            print(*row, sep='\t', file=f)


def write_read_matrix(output_dir: Path, cells: List[Cell]):
    """Create a Read-count matrix with cells as columns and clone IDs as rows"""
    clone_ids = set()
    for cell in cells:
        clone_ids.update(clone_id for clone_id in cell.counts)
    clone_ids = sorted(clone_ids)
    all_counts = [cell.counts for cell in cells]
    with open(output_dir / "read_count_matrix.csv", "w") as f:
        f.write(",")
        f.write(",".join(cell.cell_id for cell in cells))
        f.write("\n")
        for clone_id in clone_ids:
            f.write(clone_id)
            f.write(",")
            values = [lic.get(clone_id, 0) for lic in all_counts]
            f.write(",".join(str(v) for v in values))
            f.write("\n")

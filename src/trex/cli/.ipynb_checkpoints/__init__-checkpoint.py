import logging
from pathlib import Path
import shutil
from types import SimpleNamespace

from .. import __version__
from trex.utils import NiceFormatter

logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    pass


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


def add_common_arguments(parser):
    """Add arguments to an ArgumentParser"""

    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="Print some extra debugging messages",
    )

    input_group = parser.add_argument_group("Input")

    input_group.add_argument(
        "--genome-name",
        metavar="NAME",
        help="Name of the genome as indicated in 'cellranger count' run with the flag --genome. "
        "Default: Auto-detected",
        default=None,
    )
    input_group.add_argument(
        "--chromosome",
        "--chr",
        help="Name of chromosome on which cloneID is located. "
        "Default: Last chromosome in BAM file",
        default=None,
    )
    input_group.add_argument(
        "--start",
        "-s",
        help="Position of first cloneID nucleotide (1-based). Default: Auto-detected",
        type=int,
        metavar="INT",
        default=None,
    )
    input_group.add_argument(
        "--end",
        "-e",
        help="Position of last cloneID nucleotide (1-based). Default: Auto-detected",
        type=int,
        metavar="INT",
        default=None,
    )
    help = (
        "Path to Cell Ranger result directory (a subdirectory 'outs' must exist) "
        "containing sequencing of the cloneID amplicon library."
    )
    input_group.add_argument(
        "--amplicon",
        "-a",
        nargs="+",
        metavar="DIRECTORY",
        help=help + " Provide these in the same order as transcriptome datasets",
        default=None,
    )
    help = "Cell Ranger directories"
    input_group.add_argument(
        "--samples",
        help="Sample names separated by comma, in the same order as " + help,
        default=None,
    )
    input_group.add_argument(
        "--prefix",
        default=False,
        action="store_true",
        help="Add sample name as prefix to cell IDs. Default: Add as suffix",
    )

    filter_group = parser.add_argument_group("Filter settings")

    filter_group.add_argument(
        "--max-dashes",
        "-m",
        help="Maximum number of dashes that are allowed in the viral ID. Default: %(default)s",
        type=int,
        metavar="INT",
        default=10,
    )
    filter_group.add_argument(
        "--max-hamming",
        help="Maximum hamming distance allowed for two cloneIDs to be called similar. "
        "Default: %(default)s",
        type=int,
        metavar="INT",
        default=5,
    )
    filter_group.add_argument(
        "--filter-cellids",
        "-f",
        metavar="CSV",
        type=Path,
        help="CSV file containing cell IDs to keep in the analysis. "
        "This flag enables to remove cells e.g. doublets",
        default=None,
    )

    output_group = parser.add_argument_group("Output directory")

    output_group.add_argument(
        "--output",
        "-o",
        "--name",
        "-n",
        metavar="DIRECTORY",
        type=Path,
        help="Name of the run directory to be created by the program. Default: %(default)s",
        default=Path("trex_run"),
    )
    output_group.add_argument(
        "--delete",
        action="store_true",
        help="Delete the run directory if it already exists",
    )

    optional_group = parser.add_argument_group(
        "Optional output files",
        description="Use these options to enable creation "
        "of additional files in the output directory",
    )

    help = (
        "Path to the input Cell Ranger directories. "
        "There must be an 'outs' subdirectory in each of these directories."
        )
    parser.add_argument(
        "path",
        type=Path,
        nargs="+",
        metavar="DIRECTORY",
        help=help,
    )

    return SimpleNamespace(input=input_group, filter=filter_group, output=output_group)

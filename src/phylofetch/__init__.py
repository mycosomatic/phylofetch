"""phylofetch — extract phylogenetic loci from fungal genome assemblies."""

__version__ = "0.1.0"

# Public API — importable directly from `phylofetch` for use in umbrella apps.
from phylofetch.assembly_utils import get_assembly_stats, detect_assembler
from phylofetch.blast_loci_utils import extract_locus, select_best_locus_group
from phylofetch.itsx_utils import run_itsx
from phylofetch.primer_utils import find_primer_amplicons, run_primer_extraction, PRIMER_CATALOGUE
from phylofetch.project_manager import RunManager
from phylofetch.ncbi_utils import fetch_and_store, set_email
from phylofetch.busco_utils import (
    build_occupancy_matrix,
    export_sc_fastas,
    scan_busco_run,
)

__all__ = [
    "__version__",
    "get_assembly_stats",
    "detect_assembler",
    "extract_locus",
    "select_best_locus_group",
    "run_itsx",
    "find_primer_amplicons",
    "run_primer_extraction",
    "PRIMER_CATALOGUE",
    "RunManager",
    "fetch_and_store",
    "set_email",
    "build_occupancy_matrix",
    "export_sc_fastas",
    "scan_busco_run",
]

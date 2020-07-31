from os.path import join, dirname, isfile

__version__ = '0.0.7_dev'

REFBANK_SQLDB_PATH = join(dirname(__file__), 'data/refbank_genomesearch.db')
UHGG_SQLDB_PATH = join(dirname(__file__), 'data/uhgg_genomesearch.db')
PHYLOPHLAN_MARKER_PATH = join(dirname(__file__), 'data/phylophlan_marker_references.dmnd')
REFBANK_MARKER_RANKS_PATH = join(dirname(__file__), 'data/refbank_marker_ranks.txt')
UHGG_MARKER_RANKS_PATH = join(dirname(__file__), 'data/uhgg_marker_ranks.txt')
REFBANK_UNIQUE_MARKERS_PATH = join(dirname(__file__), 'data/refbank_unique_markers')
UHGG_UNIQUE_MARKERS_PATH = join(dirname(__file__), 'data/uhgg_unique_markers')

PRODIGAL_PATH = join(dirname(__file__), 'bin/prodigal.linux')
DIAMOND_PATH = join(dirname(__file__), 'bin/diamond')

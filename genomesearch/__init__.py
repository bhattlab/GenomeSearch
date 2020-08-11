from os.path import join, dirname, isfile

__version__ = '0.0.8_dev'

REFBANK_SQLDB_PATH = join(dirname(__file__), 'data/refbank_genomesearch.db')
UHGG_SQLDB_PATH = join(dirname(__file__), 'data/uhgg_genomesearch.db')
PHYLOPHLAN_MARKER_PATH = join(dirname(__file__), 'data/phylophlan_marker_references.dmnd')
REFBANK_MARKER_RANKS_ONE_PATH = join(dirname(__file__), 'data/refbank_markers_one.tsv')
REFBANK_MARKER_RANKS_ONEZERO_PATH = join(dirname(__file__), 'data/refbank_markers_onezero.tsv')
UHGG_MARKER_RANKS_ONE_PATH = join(dirname(__file__), 'data/uhgg_markers_one.tsv')
UHGG_MARKER_RANKS_ONEZERO_PATH = join(dirname(__file__), 'data/uhgg_markers_onezero.tsv')
REFBANK_UNIQUE_MARKERS_PATH = join(dirname(__file__), 'data/refbank_unique_markers')
UHGG_UNIQUE_MARKERS_PATH = join(dirname(__file__), 'data/uhgg_unique_markers')

PRODIGAL_PATH = join(dirname(__file__), 'bin/prodigal.linux')
DIAMOND_PATH = join(dirname(__file__), 'bin/diamond')
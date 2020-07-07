from os.path import join, dirname, isfile

__version__ = '0.0.1_dev'

SQLDB_PATH = join(dirname(__file__), 'data/genomesearch.db')
MARKER_RANKS_PATH = join(dirname(__file__), 'data/marker_ranks.txt')
UNIQUE_MARKERS_PATH = join(dirname(__file__), 'data/unique_markers')

PRODIGAL_PATH = join(dirname(__file__), 'bin/prodigal')
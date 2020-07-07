from genomesearch import *
import click
import wget
from os import makedirs, remove
from os.path import join, dirname, isfile
from multiprocessing import Pool
from itertools import cycle

def _download(threads, force):
    click.echo("#### INPUT PARAMETERS ####")
    try:
        num_markers = int(input("How many marker genes do you want to use? (This can be any number between 1 and 400)\n[default=150] >> ") or "150")
        if num_markers > 400 or num_markers < 1:
            raise Exception('wrong_markers')
    except:
        click.echo("ERROR!")
        click.echo("Please input a number between 1 and 400, the default is 150.")
        num_markers = int(input("How many marker genes do you want to use?]\n[default=150] >> ") or "150")
    click.echo("####################")

    if not isfile(SQLDB_PATH) or force:
        if isfile(SQLDB_PATH):
            remove(SQLDB_PATH)
        click.echo("Downloading genomesearch SQL database...")
        makedirs(dirname(SQLDB_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/genomesearch/downloads/genomesearch.db', SQLDB_PATH)
        click.echo()

    if not isfile(PHYLOPHLAN_MARKER_PATH) or force:
        if isfile(PHYLOPHLAN_MARKER_PATH):
            remove(PHYLOPHLAN_MARKER_PATH)
        click.echo("Downloading phylophlan marker gene references...")
        makedirs(dirname(PHYLOPHLAN_MARKER_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/genomesearch/downloads/phylophlan_marker_references.dmnd', PHYLOPHLAN_MARKER_PATH)
        click.echo()

    markers = []
    with open(MARKER_RANKS_PATH) as infile:
        for line in infile:
            marker = line.strip()
            markers.append(marker)

    click.echo("Downloading unique marker gene database...")
    makedirs(UNIQUE_MARKERS_PATH, exist_ok=True)
    markers = list(zip(markers[:num_markers], cycle([force])))
    with Pool(processes=threads) as pool:
        pool.starmap(download_unique_marker, markers)
    click.echo("Finished downloading...")

def download_unique_marker(marker, force):
    remote_path_dmnd = 'https://storage.googleapis.com/genomesearch/downloads/unique_markers/' + marker + '.unique.dmnd'
    remote_path_tsv = 'https://storage.googleapis.com/genomesearch/downloads/unique_markers/' + marker + '.unique.tsv'
    local_path_dmnd = join(UNIQUE_MARKERS_PATH, marker + '.unique.dmnd')
    local_path_tsv = join(UNIQUE_MARKERS_PATH, marker + '.unique.tsv')

    if not isfile(local_path_dmnd) or force:
        if isfile(local_path_dmnd):
            remove(local_path_dmnd)
        wget.download(remote_path_dmnd, local_path_dmnd)
        print()
    if not isfile(local_path_tsv) or force:
        if isfile(local_path_tsv):
            remove(local_path_tsv)
        wget.download(remote_path_tsv, local_path_tsv)
        print()

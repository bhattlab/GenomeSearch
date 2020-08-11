from genomesearch import *
import click
import wget, sys
from os import makedirs, remove
from os.path import join, dirname, isfile
from multiprocessing import Pool
from itertools import cycle


def _install(threads, force):

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

    if not isfile(REFBANK_SQLDB_PATH) or force:
        if isfile(REFBANK_SQLDB_PATH):
            remove(REFBANK_SQLDB_PATH)
        click.echo("Downloading refbank genomesearch SQL database...")
        makedirs(dirname(REFBANK_SQLDB_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/genomesearch/downloads/refbank_genomesearch.db', REFBANK_SQLDB_PATH)
        click.echo()

    if not isfile(UHGG_SQLDB_PATH) or force:
        if isfile(UHGG_SQLDB_PATH):
            remove(UHGG_SQLDB_PATH)
        click.echo("Downloading uhgg genomesearch SQL database...")
        makedirs(dirname(UHGG_SQLDB_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/genomesearch/downloads/uhgg_genomesearch.db', UHGG_SQLDB_PATH)
        click.echo()

    if not isfile(PHYLOPHLAN_MARKER_PATH) or force:
        if isfile(PHYLOPHLAN_MARKER_PATH):
            remove(PHYLOPHLAN_MARKER_PATH)
        click.echo("Downloading phylophlan marker gene references...")
        makedirs(dirname(PHYLOPHLAN_MARKER_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/genomesearch/downloads/phylophlan_marker_references.dmnd', PHYLOPHLAN_MARKER_PATH)
        click.echo()

    markers = []
    with open(REFBANK_MARKER_RANKS_ONEZERO_PATH) as infile:
        for line in infile:
            marker = line.strip().split()[0]
            markers.append(marker)

    click.echo("Downloading refbank unique marker gene database...")
    makedirs(REFBANK_UNIQUE_MARKERS_PATH, exist_ok=True)
    select_markers = list(zip(markers[:num_markers], cycle([force])))

    with Pool(processes=threads) as pool:
        pool.starmap(download_refbank_unique_marker, select_markers)
    click.echo("Finished downloading...")

    markers = []
    with open(UHGG_MARKER_RANKS_ONEZERO_PATH) as infile:
        for line in infile:
            marker = line.strip().split()[0]
            markers.append(marker)

    click.echo("Downloading uhgg unique marker gene database...")
    makedirs(UHGG_UNIQUE_MARKERS_PATH, exist_ok=True)
    if num_markers > 399:
        num_markers = 399
    select_markers = list(zip(markers[:num_markers], cycle([force])))
    with Pool(processes=threads) as pool:
        pool.starmap(download_uhgg_unique_marker, select_markers)
    click.echo("Finished installing...")


def download_refbank_unique_marker(marker, force):
    remote_path_dmnd = 'https://storage.googleapis.com/genomesearch/downloads/refbank_unique_markers/' + marker + '.unique.dmnd'
    remote_path_pkl = 'https://storage.googleapis.com/genomesearch/downloads/refbank_unique_markers/' + marker + '.unique.pkl'
    local_path_dmnd = join(REFBANK_UNIQUE_MARKERS_PATH, marker + '.unique.dmnd')
    local_path_pkl = join(REFBANK_UNIQUE_MARKERS_PATH, marker + '.unique.pkl')

    if not isfile(local_path_dmnd) or force:
        if isfile(local_path_dmnd):
            remove(local_path_dmnd)
        wget.download(remote_path_dmnd, local_path_dmnd)
        print()

    if not isfile(local_path_pkl) or force:
        if isfile(local_path_pkl):
            remove(local_path_pkl)
        wget.download(remote_path_pkl, local_path_pkl)
        print()


def download_uhgg_unique_marker(marker, force):
    print(marker, force)

    remote_path_dmnd = 'https://storage.googleapis.com/genomesearch/downloads/uhgg_unique_markers/' + marker + '.unique.dmnd'
    remote_path_pkl = 'https://storage.googleapis.com/genomesearch/downloads/uhgg_unique_markers/' + marker + '.unique.pkl'
    local_path_dmnd = join(UHGG_UNIQUE_MARKERS_PATH, marker + '.unique.dmnd')
    local_path_pkl = join(UHGG_UNIQUE_MARKERS_PATH, marker + '.unique.pkl')

    if not isfile(local_path_dmnd) or force:
        if isfile(local_path_dmnd):
            remove(local_path_dmnd)
        wget.download(remote_path_dmnd, local_path_dmnd)
        print()

    if not isfile(local_path_pkl) or force:
        if isfile(local_path_pkl):
            remove(local_path_pkl)
        wget.download(remote_path_pkl, local_path_pkl)
        print()

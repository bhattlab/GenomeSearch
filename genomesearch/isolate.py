from os import makedirs
from os.path import isdir, join
import shutil
import click
import sys
import os
from genomesearch.prodigal import run_prodigal
from genomesearch import *
from Bio import SeqIO
from subprocess import run, DEVNULL
from collections import defaultdict
import sqlite3
from glob import glob
import numpy as np
import time
from multiprocessing import Pool
import pickle

def _refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):

    tmpdir = join(outdir, 'tmp')

    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(outdir)
        makedirs(tmpdir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    prodigal_start = time.time()
    if fasta_type == 'genome':
        click.echo("Running prodigal...")
        run_prodigal(PRODIGAL_PATH, fasta, tmpdir, meta=False)
        proteome_path = join(tmpdir, 'prodigal.faa')
    elif fasta_type == 'proteome':
        proteome_path = fasta
    prodigal_end = time.time()

    marker_gene_start = time.time()
    if fasta_type == 'proteome' or fasta_type == 'genome':
        marker_output = join(outdir, prefix+'.markers.faa')
        click.echo("Identifying marker genes...")
        get_marker_genes(proteome_path, marker_output, prefix, threads)
    elif fasta_type == 'markers':
        marker_output = fasta
    marker_gene_end = time.time()

    click.echo("Searching for closest genomes in database...")
    closest_genomes_path, gene_count_time, closest_genomes_time = get_refbank_closest_genomes(
        marker_output, num_markers, tmpdir, threads, max_target_seqs, min_percent_identity
    )

    outpath = join(outdir, prefix+'.closest_genomes.tsv')
    shutil.move(closest_genomes_path, outpath)

    if not keep_intermediate:
        shutil.rmtree(tmpdir)

    print()
    print("COMPLETE.")
    print("Prodigal runtime: %f" % (prodigal_end-prodigal_start))
    print("Marker gene runtime: %f" % (marker_gene_end - marker_gene_start))
    print("Closest genome runtime: %f" % closest_genomes_time)
    print("Gene count runtime: %f" % gene_count_time)

def _uhgg(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):

    tmpdir = join(outdir, 'tmp')
    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(outdir)
        makedirs(tmpdir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    prodigal_start = time.time()
    if fasta_type == 'genome':
        click.echo("Running prodigal...")
        run_prodigal(PRODIGAL_PATH, fasta, tmpdir, meta=False)
        proteome_path = join(tmpdir, 'prodigal.faa')
    elif fasta_type == 'proteome':
        proteome_path = fasta
    prodigal_end = time.time()

    marker_gene_start = time.time()
    if fasta_type == 'proteome' or fasta_type == 'genome':
        marker_output = join(outdir, prefix+'.markers.faa')
        click.echo("Identifying marker genes...")
        get_marker_genes(proteome_path, marker_output, prefix, threads)
    elif fasta_type == 'markers':
        marker_output = fasta
    marker_gene_end = time.time()

    click.echo("Searching for closest genomes in database...")
    closest_genomes_path, gene_count_time, closest_genomes_time = get_uhgg_closest_genomes(
        marker_output, num_markers, tmpdir, threads, max_target_seqs, min_percent_identity
    )

    outpath = join(outdir, prefix+'.closest_genomes.tsv')
    shutil.move(closest_genomes_path, outpath)

    if not keep_intermediate:
        shutil.rmtree(tmpdir)

    print()
    print("COMPLETE.")
    print("Prodigal runtime: %f" % (prodigal_end-prodigal_start))
    print("Marker gene runtime: %f" % (marker_gene_end - marker_gene_start))
    print("Closest genome runtime: %f" % closest_genomes_time)
    print("Gene count runtime: %f" % gene_count_time)

def get_marker_genes(protein_fasta_path, outfile, prefix, threads):
    command = '{0} blastp --query {1} --out {2}.dmd.tsv --outfmt 6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore --db {3} --threads {4}'.format(
        DIAMOND_PATH, protein_fasta_path, outfile, PHYLOPHLAN_MARKER_PATH, threads)
    print('diamond command:', command)

    run(command.split(), stdout=DEVNULL, stderr=DEVNULL)

    top_markers = dict()
    with open(outfile + '.dmd.tsv') as infile:
        for line in infile:
            qseqid, sseqid, qlen, slen, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore =  line.strip().split('\t')
            qlen, slen, qstart, qend, sstart, send = int(qlen), int(slen), int(qstart), int(qend), int(sstart), int(send)
            length, evalue, bitscore = int(length), float(evalue), float(bitscore)
            finding = (qseqid, length, evalue, bitscore)
            marker = sseqid.split('_')[1]

            qaln = qend - qstart
            saln = send - sstart

            query_smaller = True
            if slen > qlen:
                query_smaller = False

            if finding[-2] >= 1e-6:
                continue

            if min([qlen, slen]) / max([qlen, slen]) < 0.85:
                continue

            if query_smaller:
                alnlength = qaln
            else:
                alnlength = saln

            if alnlength / min([qlen, slen]) < 0.85:
                continue

            if marker not in top_markers:
                top_markers[marker] = finding
            else:
                if finding[-1] > top_markers[marker][-1]:
                    top_markers[marker] = finding

    marker2gene = dict()
    gene2marker = dict()
    for rec in top_markers:
        marker2gene[rec] = top_markers[rec][0]
        gene2marker[top_markers[rec][0]] = rec

    records = []
    for rec in SeqIO.parse(protein_fasta_path, 'fasta'):
        if rec.id in gene2marker:
            rec.id = gene2marker[rec.id] + '__' + rec.id + '__' + prefix
            rec.description = rec.id
            records.append(rec)

    SeqIO.write(records, outfile, 'fasta')

def get_refbank_closest_genomes(marker_genes_fasta, num_markers, outdir, threads, max_target_seqs, min_percent_identity):
    closest_genomes_start = time.time()
    conn = sqlite3.connect(REFBANK_SQLDB_PATH)
    c = conn.cursor()

    c.execute("SELECT genome_id, taxon_id FROM genome;")

    genome2taxid = dict()
    for line in c.fetchall():
        genome2taxid[line[0]] = line[1]

    c.execute("SELECT taxon_id,phylum,family,species FROM taxon;")
    taxon2species = dict()
    for line in c.fetchall():
        taxon2species[line[0]] = (line[1], line[2], line[3])

    split_markers_dir = os.path.join(outdir, 'markers')
    diamond_dir = os.path.join(outdir, 'diamond')

    os.makedirs(split_markers_dir, exist_ok=True)
    os.makedirs(diamond_dir, exist_ok=True)

    marker2path = dict()
    for rec in SeqIO.parse(marker_genes_fasta, 'fasta'):
        marker = rec.id.split('__')[0]
        SeqIO.write([rec], os.path.join(split_markers_dir, marker + '.faa'), 'fasta')
        marker2path[marker] = os.path.join(split_markers_dir, marker + '.faa')

    markers = []
    with open(REFBANK_MARKER_RANKS_ONEZERO_PATH) as infile:
        count = 0
        for line in infile:
            marker = line.strip().split()[0]
            if marker in marker2path:
                markers.append(marker)
                count += 1
            if count == num_markers:
                break

    args = [(marker, split_markers_dir, diamond_dir, max_target_seqs, min_percent_identity) for marker in markers]
    with Pool(processes=threads) as pool:
        pool.starmap(run_refbank_unique_marker_search, args)

    total_markers = len(glob(diamond_dir + '/*tsv'))
    closest_genomes_end = time.time()
    gene_count_start = time.time()
    all_markers = set()
    all_pident = defaultdict(list)
    for f1 in glob(diamond_dir + '/*tsv'):

        marker = os.path.basename(f1).split('.')[0]
        all_markers.add(marker)

        seq_mapping = pickle.load(open(join(REFBANK_UNIQUE_MARKERS_PATH, marker + '.unique.pkl'), "rb"))

        with open(diamond_dir + '/' + marker + '.dmd.tsv') as infile:
            for line in infile:
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split(
                    '\t')
                pident, evalue = float(pident), float(evalue)
                if evalue >= 1e-4:
                    continue
                for genome in seq_mapping[sseqid]:
                    all_pident[genome].append((marker, pident))

    outfile = open(os.path.join(outdir, 'closest_genomes.tsv'), 'w')

    print(*(['genome', 'taxon_id', 'phylum', 'family', 'species', 'num_markers', 'total_markers', 'avg_pident'] + [marker for marker in all_markers]), sep='\t',
          file=outfile)
    closest_genomes = []
    for genome in all_pident:

        taxid = genome2taxid[int(genome)]

        if len(all_pident[genome]) / float(total_markers) < 0.25:
            continue

        marker_pident = dict(all_pident[genome])
        pidents = []
        for marker in all_markers:
            try:
                pidents.append(marker_pident[marker])
            except:
                pidents.append(None)
        closest_genomes.append(
            [genome, taxid, taxon2species[taxid][0], taxon2species[taxid][1], taxon2species[taxid][2],
             len(all_pident[genome]), total_markers, np.mean(list(marker_pident.values()))] + pidents)

    closest_genomes = list(reversed(sorted(closest_genomes, key=lambda x: x[7])))

    for res in closest_genomes:
        print(*res, sep='\t', file=outfile)

    outfile.close()

    gene_count_end = time.time()

    return os.path.join(outdir, 'closest_genomes.tsv'), gene_count_end - gene_count_start, closest_genomes_end - closest_genomes_start

def get_uhgg_closest_genomes(marker_genes_fasta, num_markers, outdir, threads, max_target_seqs, min_percent_identity):
    closest_genomes_start = time.time()
    conn = sqlite3.connect(UHGG_SQLDB_PATH)
    c = conn.cursor()

    c.execute("SELECT genome_id,genome_type,taxon_id,species_rep,completeness,contamination,country,continent  FROM genome;")

    genome2taxid = dict()
    for line in c.fetchall():
        genome2taxid[line[0]] = tuple(line[1:])

    c.execute("SELECT taxon_id,phylum,family,species FROM taxon;")
    taxon2species = dict()
    for line in c.fetchall():
        taxon2species[line[0]] = (line[1], line[2], line[3])

    split_markers_dir = os.path.join(outdir, 'markers')
    diamond_dir = os.path.join(outdir, 'diamond')

    os.makedirs(split_markers_dir, exist_ok=True)
    os.makedirs(diamond_dir, exist_ok=True)

    marker2path = dict()
    for rec in SeqIO.parse(marker_genes_fasta, 'fasta'):
        marker = rec.id.split('__')[0]
        SeqIO.write([rec], os.path.join(split_markers_dir, marker + '.faa'), 'fasta')
        marker2path[marker] = os.path.join(split_markers_dir, marker + '.faa')

    markers = []
    with open(UHGG_MARKER_RANKS_ONEZERO_PATH) as infile:
        count = 0
        for line in infile:
            marker = line.strip().split()[0]
            if marker in marker2path:
                markers.append(marker)
                count += 1
            if count == num_markers:
                break

    args = [(marker, split_markers_dir, diamond_dir, max_target_seqs, min_percent_identity) for marker in markers]
    with Pool(processes=threads) as pool:
        pool.starmap(run_uhgg_unique_marker_search, args)

    total_markers = len(glob(diamond_dir + '/*tsv'))
    closest_genomes_end = time.time()
    gene_count_start = time.time()
    all_markers = set()
    all_pident = defaultdict(list)
    for f1 in glob(diamond_dir + '/*tsv'):

        marker = os.path.basename(f1).split('.')[0]
        all_markers.add(marker)

        seq_mapping = pickle.load(open(join(UHGG_UNIQUE_MARKERS_PATH, marker + '.unique.pkl'), "rb"))

        with open(diamond_dir + '/' + marker + '.dmd.tsv') as infile:
            for line in infile:
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split(
                    '\t')
                pident, evalue = float(pident), float(evalue)
                if evalue >= 1e-4:
                    continue
                for genome in seq_mapping[sseqid]:
                    all_pident[genome].append((marker, pident))

    outfile = open(os.path.join(outdir, 'closest_genomes.tsv'), 'w')

    print(*(['genome', 'genome_type', 'species_rep', 'phylum', 'family', 'species', 'completeness', 'contamination', 'country', 'continent', 'num_markers', 'total_markers', 'avg_pident'] + [marker for marker in all_markers]), sep='\t',
          file=outfile)
    closest_genomes = []
    for genome in all_pident:
        genome_type, taxid, species_rep, completeness, contamination, country, continent = genome2taxid[genome]

        if len(all_pident[genome]) / float(total_markers) < 0.25:
            continue

        marker_pident = dict(all_pident[genome])
        pidents = []
        for marker in all_markers:
            try:
                pidents.append(marker_pident[marker])
            except:
                pidents.append(None)
        closest_genomes.append(
            [genome, genome_type, species_rep, taxon2species[taxid][0], taxon2species[taxid][1], taxon2species[taxid][2],
             completeness, contamination, country, continent,
             len(all_pident[genome]), total_markers, np.mean(list(marker_pident.values()))] + pidents)

    closest_genomes = list(reversed(sorted(closest_genomes, key=lambda x: x[12])))

    for res in closest_genomes:
        print(*res, sep='\t', file=outfile)

    outfile.close()

    gene_count_end = time.time()

    return os.path.join(outdir, 'closest_genomes.tsv'), gene_count_end - gene_count_start, closest_genomes_end - closest_genomes_start

def run_refbank_unique_marker_search(marker, split_markers_dir, diamond_dir, max_target_seqs, min_percent_identity):
    db = join(REFBANK_UNIQUE_MARKERS_PATH, marker + '.unique.dmnd')
    marker = os.path.basename(db).split('.')[0]

    command = '{0} blastp -k {1} --id {5} --query {2} --out {3}.dmd.tsv --outfmt 6 --db {4}'.format(
        DIAMOND_PATH, max_target_seqs, os.path.join(split_markers_dir, marker + '.faa'),
        os.path.join(diamond_dir, marker), db, min_percent_identity)
    print('diamond command:', command)
    run(command.split(), stdout=DEVNULL, stderr=DEVNULL)

def run_uhgg_unique_marker_search(marker, split_markers_dir, diamond_dir, max_target_seqs, min_percent_identity):
    db = join(UHGG_UNIQUE_MARKERS_PATH, marker + '.unique.dmnd')
    marker = os.path.basename(db).split('.')[0]

    command = '{0} blastp -k {1} --id {5} --query {2} --out {3}.dmd.tsv --outfmt 6 --db {4}'.format(
        DIAMOND_PATH, max_target_seqs, os.path.join(split_markers_dir, marker + '.faa'),
        os.path.join(diamond_dir, marker), db, min_percent_identity)
    print('diamond command:', command)
    run(command.split(), stdout=DEVNULL, stderr=DEVNULL)

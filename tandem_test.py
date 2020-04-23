# /usr/bin/env python

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import re
import itertools
import sys
import os
import pandas as pd
import glob
import gzip
import shutil
from ftplib import FTP

def download():

    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()

    source = ('/genomes/refseq/')
    filename = ('assembly_summary_refseq.txt')
    ftp.cwd(source)
    localfile = open(filename, 'wb')
    ftp.retrbinary("RETR " + filename, localfile.write)

    with open('species_list.txt', 'r') as list:
        with open('assembly_summary_refseq.txt', 'r') as summary:

            print('############## Downloading ################')

            list_lines = list.readlines()
            version_specie = {}

            for line in summary:
                for line2 in list_lines:
                    if (re.split('\s', line2)[0] + ' ' + re.split('\s', line2)[1]) in line and re.search(('representative genome|reference genome'), line):
                        pattern = re.compile('/genomes/all/\w{3}/\d{3}/\d{3}/\d{3}/GCF_\d+.\d+_.+\t')
                        x = pattern.findall(line)
                        vers = str(x)[2:-6]
                        version_specie.update({vers: ((re.split('[\s]', line2)[0] + '_' + re.split('[\s]', line2)[1]).strip())})

            for v, s in version_specie.items():
                print(re.split('_', s)[0] + ':', end = ' ')

                cwd = os.getcwd()
                specie = re.split('_', s)[0]
                version = re.split('[/]', v)[7]

                if not os.path.exists(cwd + '/query/' + specie + '.gtf') and not os.path.exists(cwd + '/query/' + specie + '.faa') and not os.path.exists(cwd + '/sub/' + specie + '.gtf') and not os.path.exists(cwd + '/sub/' + specie + '.faa'):
                    print('downloading...', end = ' ')
                    ftp.cwd(v)
                    gtf = (version + '_genomic.gtf.gz')
                    faa = (version + '_protein.faa.gz')
                    localfile = open(gtf, 'wb')
                    ftp.retrbinary("RETR " + gtf, localfile.write)
                    localfile = open(faa, 'wb')
                    ftp.retrbinary("RETR " + faa, localfile.write)

                    with gzip.open(version + '_genomic.gtf.gz', 'r') as gtf, open((specie + '.gtf'), 'wb') as localfile:
                        shutil.copyfileobj(gtf, localfile)
                        with gzip.open(version + '_protein.faa.gz', 'r') as faa, open((specie + '.faa'),'wb') as localfile:
                            shutil.copyfileobj(faa, localfile)
                    os.remove(version + '_genomic.gtf.gz')
                    os.remove(version + '_protein.faa.gz')
                    print('done')

                else:
                    print('already exists')

    ftp.quit()
download()

def parsing():

    def filebrowser(ext=""):
        "Returns files with an extension"
        return [f for f in glob.glob(f"*{ext}")]

    gtf = filebrowser(".gtf")
    faa = filebrowser(".faa")

    print('\n################ Parsing ##################')

    cwd = os.getcwd()

    for g in gtf:
        string = {}
        stringac = {}
        specie = os.path.splitext(g)[0]
        zprocessed = 0
        if not os.path.exists(cwd + '/query/' + specie + '_main.faa') and not os.path.exists(cwd + '/query/' + specie + '.tsv') and not os.path.exists(cwd + '/sub/' + specie + '_main.faa') and not os.path.exists(cwd + '/sub/' + specie + '.tsv'):
            with open(g) as gtf:
                for line in gtf:
                    if 'NC_' in line:
                        if 'protein_coding' in line:
                            gene = (re.split('[\t"]', line)[9])
                            chromosome = re.split('[\t]', line)[0]
                            strand = re.split('[\t]', line)[6]
                            name = (chromosome + '\t' + gene + '\t' + strand)
                        elif 'start_codon' in line:
                            if 'protein_id' in line:
                                ac = (re.search('; protein_id "(.*)', line)).group(1)
                                ac = re.split('"', ac)[0]
                                start = (re.split('\t', line))[3]
                                product = (re.search('; product "(.*)"; protein', line)).group(1)
                                stringac.setdefault(gene, []).append(ac)
                        elif 'stop_codon' in line:
                            stop = (re.split('\t', line))[4]
                            string.setdefault((name), {}).update({ac: [start, stop, product]})

            print('Chromosome\tGene\tStrand\tCDS\tStart\tStop\tProduct', file=open(specie + '.tsv', "w"))
            for e, i in string.items():
                first = list(i.items())[0]
                out = dict(itertools.islice(i.items(), 1))
                print(e + '\t' + '\t'.join("{}\t{}".format(k, '\t'.join(v)) for k, v in out.items()), file=open(specie + '.tsv', "a"))

            tab = pd.read_csv(specie + '.tsv', sep='\t', header=0)
            tab = tab.sort_values(['Chromosome', 'Start'], ascending=(True, True))
            tab.to_csv(specie + '.tsv', sep='\t', index=False)

            for f in faa:
                if os.path.splitext(g)[0] == os.path.splitext(f)[0]:
                    fasta_sequences = SeqIO.parse(open(f), 'fasta')
                    for fasta in fasta_sequences:
                        for e, i in stringac.items():
                            if i[0] in fasta.id:
                                print('>' + e + '_' + specie + '_' + fasta.id + '\n' + fasta.seq, file=open(specie + '_main.faa', 'a'))
                                zprocessed += 1
                                sys.stdout.write('\r')
                                sys.stdout.write((specie + ": {} of {} main isoforms found").format(zprocessed, len(stringac.items())))
                                sys.stdout.flush()
                    sys.stdout.write('\n')
        else:
            print(specie + ': already exists')
    else:
        print('Done')
parsing()

def movefiles():

    cwd = os.getcwd()
    dirq = os.path.join(cwd + '/query')
    dirs = os.path.join(cwd + '/sub')
    dirm = os.path.join(cwd + '/sub/main')
    if not os.path.exists(dirq):
        os.mkdir(dirq)
    if not os.path.exists(dirs):
        os.mkdir(dirs)
    if not os.path.exists(dirm):
        os.mkdir(dirm)

    queries = []
    subbies = []
    with open('species_list.txt', 'r') as list:
        list_lines = list.readlines()
        for line in list_lines:
            if 'query' in line:
                queries.append(re.split(' ', line)[0])
            if not 'query' in line:
                subbies.append(re.split(' ', line)[0])


    for line in queries:
        if not os.path.exists(cwd + '/query/' + line + '_main.faa'):
            try:
                shutil.move(cwd + '/' + line + '_main.faa', cwd + '/query/main/' + line + '_main.faa')
                shutil.move(cwd + '/' + line + '.tsv', cwd + '/query/' + line + '.tsv')
                shutil.move(cwd + '/' + line + '.faa', cwd + '/query/' + line + '.faa')
                shutil.move(cwd + '/' + line + '.gtf', cwd + '/query/' + line + '.gtf')
            except FileNotFoundError:
                pass

    for lines in subbies:
        if not os.path.exists(cwd + '/sub/' + line + '_main.faa'):
            try:
                shutil.move(cwd + '/' + lines + '_main.faa', cwd + '/sub/main/' + lines + '_main.faa')
                shutil.move(cwd + '/' + lines + '.tsv', cwd + '/sub/' + lines + '.tsv')
                shutil.move(cwd + '/' + lines + '.faa', cwd + '/sub/' + lines + '.faa')
                shutil.move(cwd + '/' + lines + '.gtf', cwd + '/sub/' + lines + '.gtf')
            except FileNotFoundError:
                pass
movefiles()

def blast_query():

    # faa = filebrowser(".faa")

    cwd = os.getcwd()
    os.chdir(cwd + '/query')

    name_seen, name2_seen = [], []
    faa_seen = {}
    print('\n############# Blast queries ###############')
    # for f in faa:
    for f in glob.glob(cwd + '/query/*.faa'): ###
        specief = re.split('/', f)[7] ###
        specie = re.split('\.', specief)[0] ###
        if not os.path.exists(cwd + '/query/' + specie + '_tandem.faa'):
            header, sequence = [], []
            fasta_sequences = SeqIO.parse(open(f), 'fasta')
            for fasta in fasta_sequences:
                header.append(fasta.id)
                sequence.append(str(fasta.seq))
            fa = dict(zip(header, sequence))

            print('Gene 1\tGene 2\tid %\tsom %\tE-value', file=open(specie + '_tandemclust' + '.tsv', 'w'))
            with open(specie + '.tsv') as tab:
                lines = [line.rstrip() for line in tab]
                res = [[lines[i], lines[i + 1]] for i in range(len(lines) - 1)]
                xprocessed, vprocessed = 0, 0
                for i in res:
                    rt = 100 * vprocessed / (len(res))
                    vprocessed += 1
                    line1, line2 = i[0], i[1]
                    [chromosome, name, strand, ac, start, stop, product] = (re.split('\t', line1))
                    [chromosome2, name2, strand2, ac2, start2, stop2, product2] = (re.split('\t', line2))
                    if strand == strand2:
                        if (strand == '+' and stop < start2) or (strand == '-' and start < stop2):
                            pair = (ac, ac2)
                            print('>' + pair[0] + '\n' + (fa.get(pair[0])), file=open(specie + "1.fa", "w"))
                            print('>' + pair[1] + '\n' + (fa.get(pair[1])), file=open(specie + "2.fa", "w"))
                            blastp_cline = blastp(query=specie + '1.fa', subject=specie + '2.fa', evalue='10e-5', max_hsps=1, out=specie + 'pair.txt')
                            stdout, stderr = blastp_cline()
                            with open(specie + 'pair.txt') as out:
                                for line in out:
                                    if 'Expect' in line:
                                        E = (re.search('Expect = (.*), ', line)).group(1)
                                    elif 'Identities' in line:
                                        Identities = (re.search('Identities = .* \((.*)\), P', line)).group(1)
                                        Positives = (re.search('Positives = .* \((.*)\), G', line)).group(1)
                                        print(name + ' (' + ac + ')' + '\t' + name2 + ' (' + ac2 + ')' + '\t' + Identities + '\t' + Positives + '\t' + E, file=open(specie + '_tandemclust' + '.tsv', 'a'))
                                        print('>' + name + '_' + specie + '_' + pair[0] + '\n' + (fa.get(pair[0])) + '\n' + '>' + name2 + '_' + specie + '_' + pair[1] + '\n' + (fa.get(pair[1])), file=open(specie + '_tandemclust' + '.faa', "a"))
                                        name_seen.append(name), name2_seen.append(name2)
                                        faa_seen.update({('>' + name + '_' + specie + '_' + pair[0] + '\n' + (fa.get(pair[0]))): ('>' + name2 + '_' + specie + '_' + pair[1] + '\n' + (fa.get(pair[1])))})
                                        xprocessed += 1
                                        sys.stdout.write('\r')
                                        sys.stdout.write((specie + ": {} pairs found [progress: {}%]").format(xprocessed, round(rt, 1)))
                                        sys.stdout.flush()
            sys.stdout.write('\n')

            with open(specie + '_tandemclust' + '.tsv', 'r') as t, open(specie + '_tandem' + '.tsv', 'w') as tn, open(specie + '_tandemclust' + '.faa', 'r') as f, open(specie + '_tandem' + '.faa', 'w') as fn:
                yprocessed = 0 - 1
                for linet in t:
                    [name, name2, id, som, Eval] = (re.split('\t', linet))
                    name_only, name2_only = re.split(' ', name)[0], re.split(' ', name2)[0]
                    if name_only not in name2_seen and name2_only not in name_seen:
                        tn.write(linet)
                        yprocessed += 1
                for k, v in faa_seen.items():
                    if (re.split('_', k)[0][1:]) in name_seen and (re.split('_', k)[0][1:]) in name2_seen:
                        continue
                    elif (re.split('_', v)[0][1:]) in name_seen and (re.split('_', v)[0][1:]) in name2_seen:
                        continue
                    else:
                        fn.write(k + '\n' + v + '\n')
                print("- " + str(yprocessed) + " isolated duplications found")

            os.remove(specie + '1' + '.fa')
            os.remove(specie + '2' + '.fa')
            os.remove(specie + 'pair' + '.txt')

        else:
            print('already exists')
    os.chdir(cwd)
blast_query()

def blast_all():

    cwd = os.getcwd()
    dir = os.path.join(cwd + '/query/tandem_vs_main')
    if not os.path.exists(dir):
        os.mkdir(dir)

    print('\n################# Blast ###################')

    sub = []
    with open('species_list.txt', 'r') as list:
        list_lines = list.readlines()
        for line in list_lines:
            if 'query' in line:
                sub.append(re.split(' ', line)[0])
    for f in glob.glob(cwd + '/query/*_tandem.faa'): ###
        specief = re.split('/', f)[7] ###
        specie = re.split('_', specief)[0] ###
        if specie in sub:
            sub.remove(specie)
        for s in sub:
            print(specie + '_tandem vs ' + s + '_main')
            cline = NcbiblastpCommandline(query=(cwd + '/query/' + specie + '_tandem.faa'), subject=(cwd + '/query/main/' + s + '_main.faa'), outfmt=7, evalue='10e-5', max_hsps=1, max_target_seqs=3, out=(cwd + '/query/tandem_vs_main/' + specie + 'tandem_' + s + '.txt'))
            cline()
blast_all()

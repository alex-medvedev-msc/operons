import numpy
import pandas
from gensim.models import Doc2Vec, Word2Vec
from gensim.models.doc2vec import TaggedDocument
from os import path, listdir
from Bio import SeqIO


def build_dict():
    aminoacids = 'ARNDCQEGHILKMFPSTWYV*XBJZ'
    aa_dict = {}
    tag = 0
    for i in range(len(aminoacids)):
        for j in range(len(aminoacids)):
            for k in range(len(aminoacids)):
                tag += 1
                aa_dict[aminoacids[i]+aminoacids[j]+aminoacids[k]] = tag
    return aa_dict


def load_proteins(data_dir):
    for file in listdir(data_dir):
        if 'proteins' in file:
            records = SeqIO.parse(path.join(data_dir, file), 'fasta')
            for record in records:
                yield (str(record.seq).rstrip('*'), record.id)


def tag_proteins(proteins):
    for i, tup in enumerate(proteins):
        protein = tup[0]
        id = tup[1]
        if i % 100000 == 0:
            print("processed {0} proteins".format(i))
        for k in range(0, 3):
            words = []
            tags = [id]
            for j in range(k, len(protein), 3):
                if len(protein) - j < 3:
                    break
                word = protein[j:j+3]
                words.append(word)
            yield TaggedDocument(words=words, tags=tags)


def load_regions(data_dir):
    for file in listdir(data_dir):
        if 'regions' in file:
            records = SeqIO.parse(path.join(data_dir, file), 'fasta')
            for record in records:
                yield (str(record.seq).rstrip('*'), record.id)


def tag_regions(regions):
    for i, tup in enumerate(regions):
        region = tup[0]
        id = tup[1]
        if i % 100000 == 0:
            print("processed {0} regions".format(i))
        for k in range(0, 9):
            words = []
            tags = [id]
            for j in range(k, len(region), 9):
                if len(region) - j < 9:
                    break
                word = region[j:j+9]
                words.append(word)
            yield TaggedDocument(words=words, tags=tags)


def proteins_from_table(filename: str):
    # gene,gi,operon_id,sequence,start,stop,strand,synonym,organism_id,genome_id
    genes = pandas.read_csv(filename)
    print(genes.shape)
    for i, gene in genes.iterrows():
        if i % 100000 == 0:
            print("processed {0} genes from table".format(i))
        for k in range(0, 3):
            words = []
            tags = [str(gene["gi"])]
            seq = gene["sequence"]
            if seq is numpy.nan:
                continue
            for j in range(k, len(seq), 3):
                if len(seq) - j < 3:
                    break
                word = seq[j:j+3]
                words.append(word)
            yield TaggedDocument(words=words, tags=tags)


def embeddings_from_fasta():
    # proteins = load_proteins('transformed')
    # documents = tag_proteins(proteins)
    # doc2vec = Doc2Vec(documents, workers=8, size=100, min_count=3, window=20, dm_concat=0)
    # doc2vec.save('models/default_tagged.doc2vec')

    regions = load_regions('transformed')
    documents = tag_regions(regions)
    doc2vec = Doc2Vec(documents, workers=8, size=50, min_count=3, window=10, dm_concat=0)
    doc2vec.save('models/default_regions_tagged.doc2vec')


def embeddings_from_table():
    documents = proteins_from_table('genes.csv')
    doc2vec = Doc2Vec(documents, workers=8, size=100, min_count=2, window=5, dm_concat=1, dbow_words=1, dm=0)
    doc2vec.save('models/doors_proteins.doc2vec')


def load_region_embeddings(regions_file, is_doc2vec=True):
    regions = pandas.read_csv(regions_file)
    regions.dropna(inplace=True)
    genes = pandas.read_csv('genes.csv')
    coding_genes = set(genes[genes.strand == '+']['gi'].values)
    K = 3
    long_count = 0
    for i, region in regions.iterrows():
        if i % 10000 == 0:
            print("processed {0} regions".format(i))
        #if region.left not in coding_genes:
        #    continue
        sequence = region.sequence
        if len(sequence) > 500:
            long_count += 1
            sequence = sequence[:250]+sequence[-250:]

        for k in range(0, K):
            words = []
            tags = ["{}_{}".format(region['left'], region['right'])]
            for j in range(k, len(sequence), K):
                if len(sequence) - j < K:
                    break
                word = sequence[j:j+K]
                words.append(word)
            if is_doc2vec:
                yield TaggedDocument(words=words, tags=tags)
            else:
                yield words
    print("count of long sequences was {}".format(long_count))


def region_embeddings():
    documents = load_region_embeddings('regions.csv')
    doc2vec = Doc2Vec(documents, workers=8, size=20, min_count=2, window=20, dm_concat=1, dbow_words=1, dm=0)
    doc2vec.save('models/doors_regions.doc2vec')


def region_word_embeddings():
    sentences = list(load_region_embeddings('regions.csv', is_doc2vec=False))
    doc2vec = Word2Vec(sentences, workers=8, size=50, min_count=2, window=10)
    doc2vec.save('models/doors_regions.word2vec')


def main():
    region_word_embeddings()


if __name__ == '__main__':
    main()

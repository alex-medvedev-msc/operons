import math
import numpy
import pandas
from gensim.models import Doc2Vec
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import seaborn
import matplotlib.pyplot as plt


def genes_to_fasta(genes, out_file):
    sequences = [SeqRecord(Seq(seq), id=str(gi), name=name, description='') for seq, gi, name in zip(genes["sequence"], genes["gi"], genes["gene"])]
    with open(out_file, "w") as file:
        SeqIO.write(sequences, file, "fasta")


def blast_correlogramm(genes_file, embeddings_file):
    embeddings = Doc2Vec.load(embeddings_file)
    genes = pandas.read_csv(genes_file)
    genes.dropna(inplace=True)
    query_sample = genes.sample(n=5000)
    subject_sample = query_sample

    query_file = "/tmp/query.fasta"
    genes_to_fasta(query_sample, query_file)
    subject_file = "/tmp/subject.fasta"
    genes_to_fasta(subject_sample, subject_file)

    out_file = "/tmp/blast_result.xml"
    # blastp = NcbiblastpCommandline(query=query_file, subject=subject_file, outfmt=5, evalue=0.001, out=out_file)
    # stdout, stderr = blastp()

    # print(stdout)
    # print(stderr)

    sims = []
    hsps = []
    with open(out_file, 'r') as result:
        results = NCBIXML.parse(result)
        for r in results:
            query_id = r.query
            for a in r.alignments:
                align_id = a.hit_id
                if query_id == align_id:
                    continue

                similarity = numpy.linalg.norm(embeddings.docvecs[query_id] - embeddings.docvecs[align_id])
                # similarity = embeddings.docvecs.similarity(query_id, align_id)
                # e = a.hsps[0].expect
                # if e < 1e-20:
                #     e = 1e-20
                # hsp = -math.log(e, 10)

                sims.append(abs(similarity))
                hsps.append(a.hsps[0].score)

    print("there are {0} hits in sample".format(len(sims)))
    print("corr coeff is {0} ".format(numpy.corrcoef(sims, hsps)))
    plt.scatter(sims, hsps)
    plt.show()


def onehot():
    genes = pandas.read_csv('genes.csv')
    data = {}
    aminoacids = {s: i for i, s in enumerate('ARNDCQEGHILKMFPSTWYV*XBJZ')}
    genes.dropna(inplace=True)
    for i, gene in genes.iterrows():
        counts = numpy.zeros((len(aminoacids),))
        for aa in gene["sequence"]:
            counts[aminoacids[aa]] += 1
        data[str(gene["gi"])] = counts/len(gene["sequence"])

    sims = []
    hsps = []
    out_file = "/tmp/blast_result.xml"
    with open(out_file, 'r') as result:
        results = NCBIXML.parse(result)
        for r in results:
            query_id = r.query
            for a in r.alignments:
                align_id = a.hit_id
                if query_id == align_id:
                    continue

                similarity = numpy.linalg.norm(data[query_id] - data[align_id])
                # similarity = embeddings.docvecs.similarity(query_id, align_id)
                # e = a.hsps[0].expect
                # if e < 1e-20:
                #     e = 1e-20
                # hsp = -math.log(e, 10)

                sims.append(abs(similarity))
                hsps.append(a.hsps[0].score)

    print("there are {0} hits in sample".format(len(sims)))
    print("corr coeff is {0} ".format(numpy.corrcoef(sims, hsps)))
    plt.scatter(sims, hsps)
    plt.show()


def simple_similarity():
    embeddings_file = 'models/doors_proteins.doc2vec'
    embeddings = Doc2Vec.load(embeddings_file)
    apaH = embeddings.docvecs.most_similar(positive=['260853263'], topn=20)
    genes = pandas.read_csv('genes.csv')

    for gene in apaH:
        print("gene: {0}, gi: {1}, similarity: {2}".format(genes[genes.gi == int(gene[0])]["gene"].values[0], gene[0], gene[1]))


def show_genes():
    genes = pandas.read_csv('genes.csv')
    #coding_strand = genes[genes['strand'] == '+']
    # gene:0,gi:1,operon_id:2,sequence:3,start:4,stop:5,strand:6,synonym:7,organism_id:8,genome_id:9
    # for i in range(1, coding_strand.values.shape[0]):
    #    if coding_strand.values[i, 4] <= coding_strand.values[i-1, 4]:
    #         print(coding_strand.values[i-1, :])
    #         print(coding_strand.values[i, :])

    print("check of coding strand finished")

    template_strand = genes[genes['strand'] == '-']
    # gene:0,gi:1,operon_id:2,sequence:3,start:4,stop:5,strand:6,synonym:7,organism_id:8,genome_id:9
    for i in range(1, template_strand.values.shape[0]):
        if template_strand.values[i, 4] <= template_strand.values[i-1, 4]:
            print(template_strand.values[i-1, :])
            print(template_strand.values[i, :])

    # print(coding_strand)


def build_dicts():
    aminoacids = 'ARNDCQEGHILKMFPSTWYV'
    aa_dict = {}
    tag = 0
    for i in range(len(aminoacids)):
        for j in range(len(aminoacids)):
            for k in range(len(aminoacids)):
                aa_dict[aminoacids[i]+aminoacids[j]+aminoacids[k]] = tag
                tag += 1
    return aa_dict, set(aminoacids)


def get_context(words: list, word_index: int, window: int):
    start = max(0, word_index - window)
    stop = min(len(words), word_index+window+1)

    context = words[start:stop]
    return context


def encode_context(context, three_grams):

    counts = numpy.zeros((len(three_grams)))
    for word in context:
        if word in three_grams:
            counts[three_grams[word]] += 1

    return counts/counts.sum()


def proteins_from_table(filename: str, limit=500):
    # gene,gi,operon_id,sequence,start,stop,strand,synonym,organism_id,genome_id
    genes = pandas.read_csv(filename)
    print(genes.shape)
    three_grams, aminoacids = build_dicts()
    encoded_contexts = []
    for i, gene in genes.sample(n=limit).iterrows():
        for k in range(0, 3):
            words = []
            seq = gene["sequence"]
            if seq is numpy.nan:
                continue
            for j in range(k, len(seq), 3):
                if len(seq) - j < 3:
                    break
                word = seq[j:j+3]
                words.append(word)

            for wi, word in enumerate(words):
                context = get_context(words, wi, 5)
                encoded = encode_context(context, three_grams)
                encoded_contexts.append(encoded)
    return numpy.array(encoded_contexts)


def analyze_pt(encoded_contexts):
    print(encoded_contexts.shape)
    means = encoded_contexts.mean(axis=0)*1000
    reversed_3grams = {tag: gram for gram, tag in build_dicts()[0].items()}
    print("most frequent 3grams in context")
    frequent_3grams = [reversed_3grams[index] for index, value in enumerate(means) if value > 0.625]
    print("number of frequent 3grams is ", len(frequent_3grams))
    print(frequent_3grams)
    print(means.shape)
    print("mean is {:4f}, std is {:4f}, min {:4f}, max {:4f}".format(means.mean(), means.std(), means.min(), means.max()))
    plt.hist(means, bins=50)
    plt.show()


def main():
    ec = proteins_from_table('genes.csv')
    analyze_pt(ec)
    # show_genes()
    # blast_correlogramm('genes.csv', 'models/doors_proteins.doc2vec')
    # onehot()

if __name__ == '__main__':
    main()

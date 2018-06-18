import math
import numpy
import pandas
from gensim.models import Doc2Vec, Word2Vec
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML
import seaborn
import pickle
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from operon import count, build_dict


def genes_to_fasta(genes, out_file):
    sequences = [SeqRecord(Seq(seq), id=str(gi), name=name, description='') for seq, gi, name in zip(genes["sequence"], genes["gi"], genes["gene"])]
    with open(out_file, "w") as file:
        SeqIO.write(sequences, file, "fasta")


def regions_to_fasta(regions, out_file):
    sequences = [SeqRecord(Seq(seq, alphabet=generic_dna), id="{}_{}".format(left, right), name='', description='') for seq, left, right
                 in zip(regions["sequence"], regions["left"], regions["right"])]
    with open(out_file, "w") as file:
        SeqIO.write(sequences, file, "fasta")


def blast_correlogramm(genes_file, embeddings_file, is_doc2vec=True):
    if is_doc2vec:
        embeddings = Doc2Vec.load(embeddings_file)
    else:
        embeddings = Word2Vec.load(embeddings_file)
    genes = pandas.read_csv(genes_file)
    genes.dropna(inplace=True)
    query_sample = genes.sample(n=2000)
    subject_sample = query_sample

    query_file = "/tmp/query.fasta"
    genes_to_fasta(query_sample, query_file)
    subject_file = "/tmp/subject.fasta"
    genes_to_fasta(subject_sample, subject_file)

    out_file = "/tmp/blast_result.xml"
    blastp = NcbiblastpCommandline(query=query_file, subject=subject_file, outfmt=5, evalue=0.001, out=out_file)
    stdout, stderr = blastp()

    # print(stdout)
    print(stderr)

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

                if is_doc2vec:
                    similarity = numpy.linalg.norm(embeddings.docvecs[query_id] - embeddings.docvecs[align_id])
                else:
                    pass
                    #for j in range(k, len(seq), 3):
                    #    if len(seq) - j < 3:
                    #        break
                    #    word = seq[j:j+3]
                    #    words.append(word)
                # similarity = embeddings.docvecs.similarity(query_id, align_id)
                e = a.hsps[0].expect
                if e < 1e-20:
                     e = 1e-20
                hsp = -math.log(e, 10)

                sims.append(abs(similarity))
                hsps.append(hsp)
                # hsps.append(a.hsps[0].score)

    print("there are {0} hits in sample".format(len(sims)))
    print("corr coeff is {0} ".format(numpy.corrcoef(sims, hsps)))

    ax = plt.subplot()
    ax.set_xlabel('Евклидово расстояние между векторами')
    ax.set_ylabel('Показатель схожести последовательностей')
    ax.scatter(sims, hsps, s=1)
    plt.show()


def build_freq_data(regions_file: str):
    regions = pandas.read_csv(regions_file)
    regions.dropna(inplace=True)
    data = {}
    combos = build_dict()

    for i, region in regions.iterrows():
        sequence = region.sequence
        # :=1 because last elem is len of sequence
        data["{}_{}".format(region["left"], region["right"])] = count(sequence, combos)[:-1]

    return data


def freq_regions(blast_file: str, region_vectors: dict):
    sims = []
    hsps = []
    with open(blast_file, 'r') as result:
        results = NCBIXML.parse(result)
        for r in results:
            query_id = r.query
            for a in r.alignments:
                align_id = a.hit_id
                if query_id == align_id:
                    continue

                similarity = numpy.linalg.norm(region_vectors[query_id] - region_vectors[align_id])
                # similarity = embeddings.docvecs.similarity(query_id, align_id)
                e = numpy.mean([hsp.expect for hsp in a.hsps])
                if e < 1e-30:
                    e = 1e-30
                hsp = -math.log(e, 10)
                hsps.append(hsp)
                sims.append(abs(similarity))
                # hsps.append(a.hsps[0].score)

    print("there are {0} hits in sample".format(len(sims)))
    print("corr coeff is {0} ".format(numpy.corrcoef(sims, hsps)))
    ax = plt.subplot()
    ax.set_xlabel('Евклидово расстояние между векторами')
    ax.set_ylabel('Показатель схожести последовательностей')
    ax.scatter(sims, hsps, s=1)
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
                e = a.hsps[0].expect
                if e < 1e-20:
                     e = 1e-20
                hsp = -math.log(e, 10)
                hsps.append(hsp)
                sims.append(abs(similarity))
                # hsps.append(a.hsps[0].score)

    print("there are {0} hits in sample".format(len(sims)))
    print("corr coeff is {0} ".format(numpy.corrcoef(sims, hsps)))
    ax = plt.subplot()
    ax.set_xlabel('Евклидово расстояние между векторами')
    ax.set_ylabel('Показатель схожести последовательностей')
    ax.scatter(sims, hsps, s=1)
    plt.show()


def blast_regions(regions_file, out_file):
    regions = pandas.read_csv(regions_file)
    regions.dropna(inplace=True)
    query_sample = regions.sample(n=5000)
    subject_sample = query_sample

    query_file = "/tmp/query.fasta"
    regions_to_fasta(query_sample, query_file)
    subject_file = "/tmp/subject.fasta"
    regions_to_fasta(subject_sample, subject_file)

    blastn = NcbiblastnCommandline(query=query_file, subject=subject_file, outfmt=5, evalue=0.1, out=out_file)
    stdout, stderr = blastn()
    print(stderr)


def blast_correlogramm_regions(embeddings_file: str, blast_file: str, region_vectors: dict, is_doc2vec=True):
    if is_doc2vec:
        embeddings = Doc2Vec.load(embeddings_file)

    sims = []
    hsps = []
    with open(blast_file, 'r') as result:
        results = NCBIXML.parse(result)
        for r in results:
            query_id = r.query
            for a in r.alignments:
                align_id = a.hit_id
                if query_id == align_id:
                    continue

                if is_doc2vec:
                    if query_id in embeddings.docvecs and align_id in embeddings.docvecs:
                        similarity = numpy.linalg.norm(embeddings.docvecs[query_id] - embeddings.docvecs[align_id])
                    else:
                        pass
                else:
                    similarity = numpy.linalg.norm(region_vectors[query_id] - region_vectors[align_id])
                    #for j in range(k, len(seq), 3):
                    #    if len(seq) - j < 3:
                    #        break
                    #    word = seq[j:j+3]
                    #    words.append(word)
                # similarity = embeddings.docvecs.similarity(query_id, align_id)
                e = numpy.mean([hsp.expect for hsp in a.hsps])
                if e < 1e-30:
                    continue
                    e = 1e-50
                hsp = -math.log(e, 10)

                sims.append(abs(similarity))
                hsps.append(hsp)
                # hsps.append(a.hsps[0].score)

    print("there are {0} hits in sample".format(len(sims)))
    print("corr coeff is {0} ".format(numpy.corrcoef(sims, hsps)))

    ax = plt.subplot()
    ax.set_xlabel('Евклидово расстояние между векторами')
    ax.set_ylabel('Показатель схожести последовательностей')
    ax.scatter(sims, hsps, s=1)
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


def build_word2vec_data(regions_file: str, embeddings_file: str):
    regions = pandas.read_csv(regions_file)
    regions.dropna(inplace=True)
    embeddings = Word2Vec.load(embeddings_file)
    K = 3
    all_vectors = {}
    table = str.maketrans(dict.fromkeys('YWRNMKSBHVD'))
    for i, region in regions.iterrows():
        if i % 10000 == 0:
            print("processed {0} regions".format(i))

        sequence = region.sequence
        if len(sequence) > 500:
            sequence = sequence[:250]+sequence[-250:]

        sequence = sequence.translate(table)

        word_vectors = []
        for k in range(0, K):
            for j in range(k, len(sequence), K):
                if len(sequence) - j < K:
                    break
                word = sequence[j:j+K]
                word_vectors.append(embeddings.wv[word])

        if len(word_vectors) == 0:
            continue

        region_matrix = numpy.array(word_vectors)
        region_vector = region_matrix.sum(axis=0)
        all_vectors["{}_{}".format(region["left"], region["right"])] = region_vector

    return all_vectors


def word2vec_freq_correlogramm(word2vec_vectors: dict, freq_vectors: dict):
    a = numpy.random.choice(list(word2vec_vectors.keys()), size=2000)
    b = numpy.random.choice(list(word2vec_vectors.keys()), size=2000)
    diffs_w2v = []
    diffs_freq = []
    for key1, key2 in zip(a, b):
        diff_w2v = numpy.linalg.norm(word2vec_vectors[key1]-word2vec_vectors[key2])
        diff_freq = numpy.linalg.norm(freq_vectors[key1]-freq_vectors[key2])
        diffs_w2v.append(diff_w2v)
        diffs_freq.append(diff_freq)

    print("corr coeff is {0} ".format(numpy.corrcoef(diffs_w2v, diffs_freq)))
    ax = plt.subplot()
    ax.set_xlabel('Евклидово расстояние между векторами word2vec')
    ax.set_ylabel('Евклидово расстояние между векторами frequency encoding')
    ax.scatter(diffs_w2v, diffs_freq, s=1)
    plt.show()


def vector_hists(vectors: dict):
    keys = numpy.random.choice(list(vectors.keys()), size=5000)

    norms = [numpy.linalg.norm(value) for key, value in vectors.items() if key in keys]

    ax = plt.subplot()
    ax.set_xlabel('Евклидово расстояние между векторами')
    ax.set_ylabel('Количество векторов')

    plt.hist(norms, bins=50)
    plt.show()


def tsne_vectors(regions_file: str, vectors: dict):
    regions = pandas.read_csv(regions_file)
    tsne = TSNE()

    in_operon = set()
    out_operon = set()
    for _, row in regions.iterrows():
        if row['in_operon']:
            in_operon.add("{}_{}".format(row["left"], row["right"]))
        else:
            out_operon.add("{}_{}".format(row["left"], row["right"]))

    keys = set(numpy.random.choice(list(vectors.keys()), size=15000))
    sample_in = [value for key, value in vectors.items() if key in keys and key in in_operon]
    sample_out = [value for key, value in vectors.items() if key in keys and key in out_operon]

    X = numpy.concatenate((numpy.array(sample_in), numpy.array(sample_out)))
    T = tsne.fit_transform(X)
    plt.scatter(T[:, 0], T[:, 1], c=["r" for _ in range(len(sample_in))] + ["b" for _ in range(len(sample_out))], s=1)
    plt.show()


def main():
    # ec = proteins_from_table('genes.csv')
    # analyze_pt(ec)
    # show_genes()
    # blast_file = '/tmp/blast_result_regions.xml'
    # blast_regions('regions.csv', blast_file)

    # word2vec_vectors = build_word2vec_data('regions.csv', 'models/doors_regions.word2vec')
    # with open('models/regions_word2vec.vectors', 'wb') as file:
        # pickle.dump(word2vec_vectors, file, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('models/regions_word2vec.vectors', 'rb') as file:
        # word2vec_vectors = pickle.load(file)
    freq_vectors = build_freq_data('regions.csv')
    tsne_vectors('regions.csv', freq_vectors)
    # vector_hists(word2vec_vectors)
    # vector_hists(freq_vectors)

    # word2vec_freq_correlogramm(word2vec_vectors, freq_vectors)

    # blast_correlogramm_regions('models/doors_regions.doc2vec', blast_file, region_vectors, False)
    # blast_correlogramm('genes.csv', 'models/doors_proteins.doc2vec')
    # onehot()
    # freq_regions(blast_file, region_vectors)

if __name__ == '__main__':
    main()

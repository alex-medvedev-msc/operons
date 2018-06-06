from gensim.models import Doc2Vec
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


def main():

    # proteins = load_proteins('transformed')
    # documents = tag_proteins(proteins)
    # doc2vec = Doc2Vec(documents, workers=8, size=100, min_count=3, window=20, dm_concat=0)
    # doc2vec.save('models/default_tagged.doc2vec')

    regions = load_regions('transformed')
    documents = tag_regions(regions)
    doc2vec = Doc2Vec(documents, workers=8, size=50, min_count=3, window=10, dm_concat=0)
    doc2vec.save('models/default_regions_tagged.doc2vec')


if __name__ == '__main__':
    main()

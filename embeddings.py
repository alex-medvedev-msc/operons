#from gensim.models import word2vec
import pandas
import numpy
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import pairwise_distances
import matplotlib.pyplot as plt


def main():

    def join(s):
        return " ".join(s)

    path = "all_gis.csv"
    operons = pandas.read_csv(path)

    raw_sentences = operons[operons.gene != '-'].groupby(by="operon_id")["gene"].aggregate(join)
    sentences = [s.split(' ') for s in raw_sentences.values]
    model = word2vec.Word2Vec(sentences, size=50, window=10, min_count=1, workers=4)
    model.save("escherichia_embeddings")
    #print(model.wv.similarity('prpC', 'lacZ'))


def do_vectors(path, embeddings_path):

    operons = pandas.read_csv(path)
    model = word2vec.Word2Vec.load(embeddings_path)

    def vectorize(operon):
        array = numpy.vstack([model.wv[gene] for gene in operon])
        return array.T.mean(axis=1)

    grouped = operons[operons.gene != '-'].groupby(by="operon_id")["gene"]
    vectorized = numpy.vstack([vectorize(group).T for name, group in grouped])
    df = pandas.DataFrame(data=vectorized, index=[name for name, _ in grouped], columns=["f{0}".format(i) for i in range(50)])
    df.to_csv("vectorized_gis.csv", float_format='%.7f', index=True, index_label='operon_id')
    #numpy.save("vectorized_gis", vectorized)


def cluster(vectors_path):
    X = pandas.read_csv(vectors_path, index_col='operon_id')
    scan = DBSCAN(min_samples=3, metric='cosine', n_jobs=2, eps=0.5)
    labels = scan.fit_predict(X.values)
    labels_frame = pandas.DataFrame(data=labels, index=X.index, columns=["labels"])
    labels_frame.to_csv("clusters_gis.csv")


def hist(vectors_path):
    X = pandas.read_csv(vectors_path, index_col='operon_id')
    idx = numpy.random.choice(X.shape[0], size=(1000,), replace=False)

    sample = X.values[idx, :]
    distances = pairwise_distances(sample, metric="cosine").flatten()
    plt.hist(distances, bins=100)
    plt.show()


if __name__ == '__main__':
    #do_vectors("all_gis.csv", "escherichia_embeddings")
    cluster("vectorized_gis.csv")

    #hist("vectorized_gis.csv")
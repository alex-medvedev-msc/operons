from gensim.models import word2vec
import pandas
import numpy
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import pairwise_distances
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize


def clear_data(path):
    operons = pandas.read_csv(path)
    grouped = operons.groupby(by="operon_id")
    clean = []
    for name, group in grouped:
        if "-" not in group["gene"].values and len(group["gene"].values) > 1:
            clean.append(name)
    clean_operons = operons[operons["operon_id"].isin(clean)]
    clean_operons.to_csv("all_gis_clean.csv")


def embeddings():

    def join(s):
        return " ".join(s)

    path = "all_gis_clean.csv"
    operons = pandas.read_csv(path)

    raw_sentences = operons.groupby(by="operon_id")["gene"].aggregate(join)
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
    scan = DBSCAN(min_samples=5, metric='l2', n_jobs=4, eps=0.3, algorithm='kd_tree')
    labels = scan.fit_predict(normalize(X.values, "l2"))
    labels_frame = pandas.DataFrame(data=labels, index=X.index, columns=["labels"])
    labels_frame.to_csv("clusters_gis.csv")


def hist(vectors_path):
    X = pandas.read_csv(vectors_path, index_col='operon_id')
    idx = numpy.random.choice(X.shape[0], size=(2000,), replace=False)

    sample = X.values[idx, :]
    distances = pairwise_distances(normalize(sample, "l2"), metric="euclidean").flatten()
    plt.hist(distances, bins=100)
    plt.show()


def plot_clusters():
    labels = pandas.read_csv("clusters_gis.csv")
    clusters = {}
    for index, row in labels[labels["labels"] >= 0].iterrows():
        if str(row["labels"]) in clusters:

            clusters[str(row["labels"])].append(str(row["operon_id"]))
        else:
            clusters[str(row["labels"])] = [str(row["operon_id"])]

    import json
    #with open("clusters.json", "w") as file:
    #    json.dump(clusters, file)

    genes = pandas.read_csv("all_gis_clean.csv").groupby(by="operon_id")
    result = {}
    for key, cluster in clusters.items():
        result[key] = {}
        for operon in cluster:
            result[key][operon] = genes.get_group(int(operon))["gene"].tolist()

    with open("operons_clusters.json", "w") as file:
        json.dump(result, file)


    #counts = labels.groupby(by="labels").size().reset_index(name='counts')
    #plt.hist(counts[counts["labels"] > 0].counts.values)
    #plt.show()


if __name__ == '__main__':
    #clear_data("all_gis.csv")
    #embeddings()
    #do_vectors("all_gis_clean.csv", "escherichia_embeddings")
    #hist("vectorized_gis.csv")

    cluster("vectorized_gis.csv")
    plot_clusters()

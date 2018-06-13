import numpy
import pandas
from gensim.models import Doc2Vec, Word2Vec
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import log_loss, accuracy_score
from sklearn.metrics import classification_report


def build_data(genes_file, embeddings_file):

    genes = pandas.read_csv(genes_file)
    embeddings = Doc2Vec.load(embeddings_file)

    strands = genes.groupby(by=["organism_id", "strand"])
    features = []
    Y = []
    missing_count = 0
    print("there are {0} strands".format(len(strands)))
    # gene:0,gi:1,operon_id:2,sequence:3,start:4,stop:5,strand:6,synonym:7,organism_id:8,genome_id:9
    for key, strand in strands:
        values = strand.values
        for i in range(1, values.shape[0]):
            gi1 = str(values[i-1, 1])
            gi2 = str(values[i, 1])
            if gi1 not in embeddings.docvecs:
                missing_count += 1
                continue
            if gi2 not in embeddings.docvecs:
                continue
            vec1 = embeddings.docvecs[gi1]
            vec2 = embeddings.docvecs[gi2]

            feature_vec = numpy.concatenate((vec1, vec2), axis=0)
            features.append(feature_vec.T)

            if values[i-1, 2] == values[i, 2]:
                # they are from the same operon
                Y.append(1)
            else:
                # they are from different operons
                Y.append(0)
        print("key {0} finished, processed {1} rows".format(key, values.shape[0]))

    X = numpy.array(features)
    print("X has shape {0}".format(X.shape))
    print("there are {0} genes missing from doc2vec".format(missing_count))
    data = pandas.DataFrame(data=X, columns=["left{0}".format(i) for i in range(100)] + ["right{0}".format(i) for i in range(100)])
    data.loc[:, "Y"] = Y
    data.to_csv("data.csv", index=False)


def count(sequence, combos):
    if len(sequence) > 500:
        sequence = sequence[:250]+sequence[-250:]

    table = str.maketrans(dict.fromkeys('YWRNMKSBHVD'))
    sequence = sequence.translate(table)
    counts = numpy.zeros((4,))
    for i in range(len(sequence)):
        counts[combos[sequence[i]]] += 1

    counts2 = numpy.zeros((16,))
    for i in range(1, len(sequence)-1):
        counts2[combos[sequence[i-1:i+1]]-4] += 1

    counts3 = numpy.zeros((64,))
    for i in range(2, len(sequence)-1):
        counts3[combos[sequence[i-2:i+1]]-4-16] += 1

    counts = counts/len(sequence)
    if sum(counts2) != 0:
        counts2 = counts2/sum(counts2)
    if sum(counts3) != 0:
        counts3 = counts3/sum(counts3)
    return numpy.concatenate((counts, counts2, counts3, numpy.array([len(sequence)])))


def build_dict():
    nc = 'ACGT'

    combos = {n: i for i, n in enumerate(nc)}
    for i in range(4):
        combos[nc[i]] = i
        for j in range(4):
            combos[nc[i]+nc[j]] = 4+i*4+j
            for k in range(4):
                combos[nc[i]+nc[j]+nc[k]] = 20+i*16+j*4+k

    return combos


def build_onehot_data(regions_file: str, genes_file: str, out_file: str):
    regions = pandas.read_csv(regions_file)
    regions.dropna(inplace=True)

    genes = pandas.read_csv(genes_file)
    # coding_genes = set(genes[genes.strand == '+']['gi'].values)

    print("regions shape is {}".format(regions.shape))

    combos = build_dict()

    features = []
    Y = []
    for i, region in regions.iterrows():
        # if region["left"] not in coding_genes:
        #    continue
        region_id = "{}_{}".format(region["left"], region["right"])

        feature_vec = count(region.sequence, combos)
        # feature_vec = embeddings.docvecs[region_id]
        features.append(feature_vec)
        Y.append(region["in_operon"])

    X = numpy.array(features)
    print("X has shape {0}".format(X.shape))
    data = pandas.DataFrame(data=X, columns=["c{0}".format(i) for i in range(4+16+64)] + ["distance"])
    data.loc[:, "Y"] = Y
    data.to_csv(out_file, index=False)


def build_regions_data(regions_file, embeddings_file):

    regions = pandas.read_csv(regions_file)
    embeddings = Doc2Vec.load(embeddings_file)
    regions.dropna(inplace=True)

    genes = pandas.read_csv('genes.csv')
    coding_genes = set(genes[genes.strand == '+']['gi'].values)

    print("regions shape is {}".format(regions.shape))

    features = []
    Y = []
    # in_operon:0,left:1,right:2,sequence:3,organism_id:4,genome_id:5

    for i, region in regions.iterrows():
        if region["left"] not in coding_genes:
            continue
        region_id = "{}_{}".format(region["left"], region["right"])
        feature_vec = embeddings.docvecs[region_id]
        features.append(feature_vec)
        Y.append(region["in_operon"])

    X = numpy.array(features)
    print("X has shape {0}".format(X.shape))
    data = pandas.DataFrame(data=X, columns=["c{0}".format(i) for i in range(20)])
    data.loc[:, "Y"] = Y
    data.to_csv("regions_data.csv", index=False)


def classify(data_file):
    data = pandas.read_csv(data_file)
    X = data.values[:, :200]
    Y = data.values[:, 200]
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, train_size=0.7)
    clf = GradientBoostingClassifier(subsample=0.7, max_depth=5, verbose=1, n_estimators=50)
    clf.fit(X_train, Y_train)
    Y_pred = clf.predict_proba(X_test)
    Y_pred_accurate = clf.predict(X_test)
    loss, accuracy = log_loss(Y_test, Y_pred), accuracy_score(Y_test, Y_pred_accurate)
    print("loss is {0:4f}, accuracy is {1:4f}".format(loss, accuracy))
    report = classification_report(Y_test, Y_pred_accurate, target_names=["out of operon", "in operon"])
    print(report)


def classify_regions(data_file: str, validate_file: str):
    data = pandas.read_csv(data_file)
    X = data.values[:, :data.shape[1]-1]
    Y = data.values[:, -1]
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, train_size=0.7)
    clf = GradientBoostingClassifier(learning_rate=0.1, subsample=0.7, max_depth=5, verbose=1, n_estimators=60)
    clf.fit(X_train, Y_train)
    Y_pred = clf.predict_proba(X_test)
    Y_pred_accurate = clf.predict(X_test)
    ones = sum(Y_pred_accurate)/len(Y_pred_accurate)
    zeros = 1 - ones
    print("number of zeros {}, ones {}".format(zeros, ones))
    loss, accuracy = log_loss(Y_test, Y_pred), accuracy_score(Y_test, Y_pred_accurate)
    print("loss is {0:4f}, accuracy is {1:4f}".format(loss, accuracy))
    report = classification_report(Y_test, Y_pred_accurate, target_names=["out of operon", "in operon"])
    print(report)
    print(clf.feature_importances_)
    # combos = {value: key for key, value in build_dict().items()}
    # for i, fi in enumerate(clf.feature_importances_[:-1]):
    #     print("for motif {} importance is {:4f}".format(combos[i], fi))

    if validate_file == "":
        return
    val_data = pandas.read_csv(validate_file)
    X_val = val_data.values[:, :val_data.shape[1]-1]
    Y_val = val_data.values[:, -1]
    Y_val_accurate = clf.predict(X_val)
    val_accuracy = accuracy_score(Y_val, Y_val_accurate)
    print("validation accuracy is {}".format(val_accuracy))
    report = classification_report(Y_val, Y_val_accurate, target_names=["out of operon", "in operon"])

    print(report)


def build_word2vec_data(regions_file: str, embeddings_file: str):
    regions = pandas.read_csv(regions_file)
    regions.dropna(inplace=True)
    embeddings = Word2Vec.load(embeddings_file)
    K = 3
    all_vectors = []
    Y = []
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

        all_vectors.append(region_vector)
        Y.append(region["in_operon"])

    X = numpy.array(all_vectors)
    print("X has shape {0}".format(X.shape))
    data = pandas.DataFrame(data=X, columns=["c{0}".format(i) for i in range(X.shape[1])])
    data.loc[:, "Y"] = Y
    data.to_csv("regions_word2vec_data.csv", index=False)
    return all_vectors


def main():
    # build_data('genes.csv', 'models/doors_proteins.doc2vec')
    # classify('data.csv')
    # out_file = 'regions_bacillus_data.csv'
    # build_onehot_data('regions_bacillus.csv', 'genes_bacillus.csv', out_file)
    # build_regions_data('regions.csv', 'models/doors_regions.doc2vec')
    # classify_regions('regions_data_oh3.csv', out_file)

    # build_word2vec_data('regions.csv', 'models/doors_regions.word2vec')

    classify_regions('regions_word2vec_data.csv', "")

if __name__ == '__main__':
    main()

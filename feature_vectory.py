import pandas
import numpy

def operon_len(operon):
    return operon.end.max() - operon.start.min()


def operon_intergenetic_mean(operon):

    if len(operon) <= 1:
        return 0

    total = 0.0
    for i in range(len(operon) - 1):
        # start = 0, end = 1
        total += operon.iat[i+1, 0] - operon.iat[i, 1]

    return total/(len(operon) - 1)


def operon_intergenetic_max(operon):

    if len(operon) <= 1:
        return 0

    max = 0
    for i in range(len(operon) - 1):
        # start = 0, end = 1
        distance = operon.iat[i+1, 0] - operon.iat[i, 1]
        if distance > max:
            max = distance

    return max


def operon_pairs_count(operon):

    if len(operon) <= 1:
        return pandas.DataFrame()

    pairs = {}
    for i in range(len(operon) - 1):
        left, right = operon.iat[i, 2], operon.iat[i+1, 2]
        if not left.startswith("lac") or not right.startswith("lac"):
            continue
        strand = operon.iat[i, 4]
        if strand == '-':
            pair = right[:4] + "_" + left[:4]
        else:
            pair = left[:4] + "_" + right[:4]
        if pair not in pairs:
            pairs[pair] = 1
        else:
            pairs[pair] += 1

    return pairs


def vector(operon):
    sorted = operon.sort_values(by=["start"])

    pairs = operon_pairs_count(sorted)
    pairs["length"] = operon_len(sorted)
    pairs["mean_distance"] = operon_intergenetic_mean(sorted)
    pairs["max_distance"] = operon_intergenetic_max(sorted)
    return pandas.DataFrame(pairs, index=sorted.operon_id)


def process_table(table_path):
    positive_path = table_path
    positive = pandas.read_csv(positive_path).loc[:, ["start", "end", "gene", "operon_id", "strand"]]
    print(len(numpy.unique(positive.operon_id.values)))
    positive_features = positive.groupby('operon_id').apply(vector).groupby('operon_id').agg(lambda x: x.iloc[0])
    return positive_features


def main():
    positive_features = process_table("table.csv")
    negative_features = process_table("table_negative.csv")
    positive_cols = positive_features.columns.values
    negative_features = negative_features.loc[:, positive_cols]
    all = positive_features.append(negative_features, ignore_index=True)
    all.loc[:len(positive_features), "y"] = 1
    all.loc[len(positive_features):, "y"] = 0
    all.fillna(0, inplace=True)
    all.to_csv("features.csv", index=False)


if __name__ == '__main__':
    main()
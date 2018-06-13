import numpy
import pandas
import json
import re


def parse_organisms(org_path, pattern):
    import json
    ids = []
    with open(org_path, 'r') as file:
        data = json.loads(file.read())

    plasmid_count = 0
    chromosome_count = 0
    for organism in data['aaData']:
        if pattern not in organism[0]:
            continue
        matches = re.findall('id=[0-9]+', organism[1])
        names = re.findall('nc=\w+', organism[1])
        types = re.findall('\(.\)', organism[1])
        for group, name, tp in zip(matches, names, types):
            ids.append((group[3:], name[3:]))
            if tp == '(P)':
                plasmid_count+=1
            else:
                chromosome_count+=1

    print('there are {} plasmids and {} chromosomes'.format(plasmid_count, chromosome_count))

    return ids


def main():
    genes = pandas.read_csv('genes.csv')
    print('number of unique genomes is {}'.format(len(numpy.unique(genes.organism_id.values))))
    grouped = genes.loc[:, ['genome_id', 'operon_id', 'gi']].groupby(by='genome_id').agg(lambda x: len(numpy.unique(x)))

    print('count of unannotated genes {}'.format(genes[genes.gene == '-'].count()))

    counts = genes.loc[:, ['operon_id', 'gi']].groupby(by='operon_id').agg('count')
    print('operons count with 2+ genes: {}, 1 gene: {}'.format(len(counts.gi[counts.gi > 1]), sum(counts.gi[counts.gi == 1])))

    print('average number of operons in chromosomes {}'.format(grouped[grouped.operon_id > 1000].operon_id.mean()))
    print('average number of gi in chromosomes {}'.format(grouped[grouped.operon_id > 1000].gi.mean()))
    print('average number of operons in plasmids {}'.format(grouped[grouped.operon_id < 1000].operon_id.mean()))
    print('average number of gi in plasmids {}'.format(grouped[grouped.operon_id < 1000].gi.mean()))
    # print('unique counts: {}'.format(grouped))
    parse_organisms('../organisms.json', 'Escherichia')

if __name__ == '__main__':
    main()
import requests
import pandas
import re
from time import time, sleep
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_link_value(link):
    match = re.search(">[0-9]+<", link)
    if match:
        return match.group(0)[1:-1]
    else:
        return ""


def load_sequences(genome_id):

    handle = Entrez.esearch("nuccore", genome_id)
    records = Entrez.read(handle)
    ids = records["IdList"]
    handle.close()

    genome_handle = Entrez.efetch("nuccore", id=ids[0],  rettype="fasta", retmode="text")
    records = list(SeqIO.parse(genome_handle, format='fasta'))
    genome_handle.close()
    if len(records) > 1:
        print("{0} in {1} ".format(len(records), [record.id for record in records]))
    return records[0]


def load_organism(id, genome_id):
    url = "http://csbl.bmb.uga.edu/DOOR/NCoperon_ajax.php?mode=DataTable&id={0}&_=1516688758122".format(id)
    genes = []

    try:
        genome = load_sequences(genome_id)
        response = requests.get(url, timeout=30)
        data = response.json()["aaData"]
        for operon in data:
            operon_id = get_link_value(operon[1])
            gis = operon[2].split("<br>")
            gene_names = operon[4].split("<br>")
            synonyms = operon[3].split("<br>")
            starts = operon[5].split("<br>")
            stops = operon[6].split("<br>")
            strand = operon[7].split("<br>")[0]
            for gi, name, synonym, start, stop in zip(gis, gene_names, synonyms, starts, stops):
                g = {}
                g["operon_id"] = operon_id
                g["gi"] = get_link_value(gi)
                g["synonym"] = synonym
                g["gene"] = name
                g["start"] = int(start)
                g["stop"] = int(stop)
                g["strand"] = strand
                if strand == "+":
                    g["sequence"] = str(genome.seq[g["start"]-1:g["stop"]+1].translate()).strip("*")
                else:
                    s3 = genome.seq[g["start"]:g["stop"]].reverse_complement()
                    g["sequence"] = str(s3.translate()).strip("*")

                genes.append(g)
    except Exception as e:
        print(e)
        return None

    return genes


def load_organism_regions(genome_id, genes: pandas.DataFrame):
    regions = []
    coding_strand = genes[(genes.genome_id == genome_id) & (genes.strand == '+')]
    template_strand = genes[(genes.genome_id == genome_id) & (genes.strand == '-')]

    genome = load_sequences(genome_id)

    # gene:0,gi:1,operon_id:2,sequence:3,start:4,stop:5,strand:6,synonym:7,organism_id:8,genome_id:9
    values = coding_strand.values
    for i in range(1, values.shape[0]):
        gi1 = str(values[i-1, 1])
        gi2 = str(values[i, 1])
        end = values[i-1, 5]
        start = values[i, 4]
        in_operon = int(values[i-1, 2] == values[i, 2])
        seq = str(genome.seq[end:start])
        region = {"sequence": seq, "left": gi1, "right": gi2, "in_operon": in_operon}
        regions.append(region)

    values = template_strand.values
    for i in range(1, values.shape[0]):
        gi1 = str(values[i-1, 1])
        gi2 = str(values[i, 1])
        end = values[i-1, 5]
        start = values[i, 4]

        in_operon = int(values[i-1, 2] == values[i, 2])
        seq = str(genome.seq[end:start].reverse_complement())
        region = {"sequence": seq, "left": gi2, "right": gi1, "in_operon": in_operon}

        regions.append(region)

    return regions


def parse_organisms(org_path, pattern, limit=0):
    import json
    ids = []
    with open(org_path, 'r') as file:
        data = json.loads(file.read())

    count = 0
    if limit == 0:
        limit = 100000
    for organism in data['aaData']:
        if pattern not in organism[0]:
            continue
        if count >= limit:
            continue
        count += 1
        matches = re.findall('id=[0-9]+', organism[1])
        names = re.findall('nc=\w+', organism[1])
        for group, name in zip(matches, names):
            ids.append((group[3:], name[3:]))

    return ids


def load_proteins(out_file: str, pattern: str, limit: int=0):
    org_path = "../organisms.json"
    ids = parse_organisms(org_path, pattern=pattern, limit=limit)

    table = pandas.DataFrame()
    for id in ids:
        now = time()
        organism_id, genome_id = id[0], id[1]
        for attempt in range(3):
            gene_entries = load_organism(organism_id, genome_id)
            if gene_entries is not None:
                break
            print("attempt {0} was not successful".format(attempt))
            sleep((attempt+1)*90)

        if gene_entries is None:
            print("organism {0} was not loaded properly, DOOR is shit".format(id))
            sleep(300)
            continue

        temp = pandas.DataFrame(gene_entries)
        temp.loc[:, "organism_id"] = organism_id
        temp.loc[:, "genome_id"] = genome_id
        table = table.append(temp, ignore_index=True)
        print("organism {0} loaded {1} gis in {2:4f}s".format(id, len(gene_entries), time() - now))
        sleep(90)

    table.to_csv(out_file, index=False)


def load_regions(genes_file: str, out_file: str, pattern: str, limit: int=0):
    org_path = "../organisms.json"
    ids = parse_organisms(org_path, pattern=pattern, limit=limit)
    genes = pandas.read_csv(genes_file)
    table = pandas.DataFrame()
    for id in ids:
        now = time()
        organism_id, genome_id = id[0], id[1]
        regions = load_organism_regions(genome_id, genes)
        if len(regions) == 0:
            continue
        temp = pandas.DataFrame(regions)
        temp.loc[:, "organism_id"] = organism_id
        temp.loc[:, "genome_id"] = genome_id
        table = table.append(temp, ignore_index=True)
        print("organism {0} loaded {1} intergenetic regions in {2:4f}s".format(id, len(regions), time() - now))

    table.to_csv(out_file, index=False)


def main():
    Entrez.email = "onl1ner.as@gmail.com"
    load_proteins("genes_bacillus.csv", "Bacillus", limit=10)
    load_regions("genes_bacillus.csv", "regions_bacillus.csv", "Bacillus coagulans 36D1", limit=10)

if __name__ == '__main__':
    main()

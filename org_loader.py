import requests
import pandas
import re
from time import time


def get_link_value(link):
    match = re.search(">[0-9]+<", link)
    if match:
        return match.group(0)[1:-1]
    else:
        return ""


def load_organism(id):
    url = "http://csbl.bmb.uga.edu/DOOR/NCoperon_ajax.php?mode=DataTable&id={0}&_=1516688758122".format(id)
    genes = []

    try:
        response = requests.get(url)
        data = response.json()["aaData"]
        for operon in data:
            operon_id = get_link_value(operon[1])
            gis = operon[2].split("<br>")
            gene_names = operon[4].split("<br>")
            for gi, name in zip(gis, gene_names):
                g = {}
                g["operon_id"] = operon_id
                g["gi"] = get_link_value(gi)
                g["gene"] = name
                genes.append(g)
    except Exception as e:
        print(e)

    return genes


def parse_organisms(path, pattern):
    import json
    org_path = "organisms.json"
    ids = []
    with open(org_path, 'r') as file:
        data = json.loads(file.read())

    for organism in data['aaData']:
        if pattern not in organism[0]:
            continue
        matches = re.findall('id=[0-9]+', organism[1])
        for group in matches:
                ids.append(group[3:])

    return ids


def main():

    org_path = "organisms.json"
    ids = parse_organisms(org_path, pattern="")[500:1000]

    table = pandas.DataFrame()
    for id in ids:
        now = time()
        gene_entries = load_organism(id)
        temp = pandas.DataFrame(gene_entries)
        temp.loc[:, "organism_id"] = id
        table = table.append(temp, ignore_index=True)
        print("organism {0} loaded {1} gis in {2}".format(id, len(gene_entries), time() - now))

    table.to_csv("all_gis_500_1000.csv", index=False)


if __name__ == '__main__':
    main()

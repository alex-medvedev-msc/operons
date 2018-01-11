import requests
import pandas
from bs4 import BeautifulSoup


def load_operon(id):
    url = "http://csbl.bmb.uga.edu/DOOR/operon.php?id={0}#basic".format(id)
    response = requests.get(url)
    html = BeautifulSoup(response.text)
    table = html.find(id="basic").select(".small_data").tbody
    operons = []
    for tr in table.find_all("tbody"):
        operon = {}
        operon["gi"] = tr.descendants[0].a.data
        operon["start"] = tr.descendants[1].data
        operon["end"] = tr.descendants[2].data
        operon["strand"] = tr.descendants[3].data
        operon["gene"] = tr.descendants[4].data
        operon["description"] = tr.descendants[7]
        operons.append(operon)
    return operons

def load_gi(gi):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&rettype=fasta&retmode=text".format(gi)

def main():
    ids_path = "ids.txt"
    ids = []
    with open(ids_path, 'r') as file:
        ids = file.readlines()

    load_operon()

    table = pandas.DataFrame()
    for i, id in enumerate(ids):
        gene_entry = load_operon(id)
        table.loc[i] = gene_entry

    table.to_csv(index=False)


if __name__ == '__main__':
    main()
import requests
import pandas
from bs4 import BeautifulSoup


def load_operon(id):
    url = "http://csbl.bmb.uga.edu/DOOR/operon.php?id={0}#basic".format(id)
    html = ""
    try:
        response = requests.get(url)
        html = BeautifulSoup(response.text, "lxml")
        table = html.find(id="basic").select(".small_data")[0]
        genes = []
        for tr in table.tbody.find_all("tr"):
            gene = {}
            tds = tr.find_all("td")
            gene["gi"] = tds[0].text
            gene["start"] = tds[1].text
            gene["end"] = tds[2].text
            gene["strand"] = tds[3].text
            gene["gene"] = tds[4].text
            gene["description"] = tds[7].text
            genes.append(gene)
        return genes
    except Exception as e:
        print(e)
        print(response.status_code)
        print(html)


def load_gi(gi):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&rettype=fasta&retmode=text".format(gi)
    fasta = requests.get(url).text
    return "".join(fasta.splitlines(False)[1:])


def main():

    ids_path = "negative.csv"
    ids = []
    with open(ids_path, 'r') as file:
        ids = [line.strip() for line in file.readlines()]

    table = pandas.DataFrame()
    index = 0
    for id in ids:
        gene_entries = load_operon(id)
        for entry in gene_entries:
            entry["sequence"] = load_gi(entry["gi"])
            entry["operon_id"] = id
            table = table.append(entry, ignore_index=True)
        print("operon {0} loaded {1} gis".format(id, len(gene_entries)))

    table.to_csv("table_negative.csv", index=False)


if __name__ == '__main__':
    main()
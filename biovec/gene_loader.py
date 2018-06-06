import pandas
import requests
import re
from time import sleep


def load_gis(filename: str):
    gis = pandas.read_csv(filename)
    return gis


def parse_gene(html):
    matches = re.findall('<pre>.*<pre>', html, re.DOTALL)
    for group in matches:
        return group.replace('<pre>', '').replace('\n', '')


def load_gene(gi):
    """http://csbl.bmb.uga.edu/DOOR/genedetail.php?id=260763803"""
    url = 'http://csbl.bmb.uga.edu/DOOR/genedetail.php?id={0}'.format(gi)
    try:
        response = requests.get(url)
        return parse_gene(response.text)
    except Exception as e:
        print(e)


def load_door_genes(gis: pandas.DataFrame):
    count = 0
    for i, row in gis.iterrows():

        if 'sequence' in row:
            continue

        gi = row["gi"]
        sequence = load_gene(gi)
        if sequence is None:
            break

        count += 1
        if count % 10 == 0:
            print('loaded {0} genes'.format(count))

        sleep(0.5)
        gis.loc[i, "sequence"] = sequence

    return gis


def main():
    gis = load_gis('../all_gis.csv')
    loaded = load_door_genes(gis)
    loaded.to_csv(index=False, path_or_buf='gis.csv')


if __name__ == '__main__':
    main()

import gzip
from Bio import SeqIO, SeqRecord
import pandas

"""
feature	class	        assembly	    assembly_unit seq_type  chromosome	genomic_accession start	end	strand	product_accession	non-redundant_refseq	related_accession	name	symbol	GeneID	locus_tag	feature_interval_length	product_length	attributes
gene	protein_coding	GCF_000006665.1	Primary       Assembly	chromosome	NC_002655.2	      190	273	+		Z_RS00005	        84		old_locus_tag=Z0001
"""


def save_proteins(filename, proteins):
    SeqIO.write()


def translate_record(source: SeqRecord, start, end, strand):
    if strand == "+":
        dna = source.seq[start:end].translate()
    else:
        dna = source[start:end].reverse_complement().translate()
    return dna


def translate(records, table: pandas.DataFrame):
    chromosome = records[0]
    plasmid = records[1]
    if len(records) > 2:
        raise Exception("records have len({0}) instead of 2".format(len(records)))

    proteins = []

    for _, gene in table.iterrows():
        if gene["chromosome"] == "chromosome":
            protein = translate_record(chromosome, gene["start"], gene["end"], gene["strand"])
        else:
            protein = translate_record(plasmid, gene["start"], gene["end"], gene["strand"])

        proteins.append(protein)

    return proteins

def load_feature_table(filename):
    with gzip.open(filename, mode='rt') as file:
        table = pandas.read_csv(file)
        return table


def ungzip(filename):
    records = []
    with gzip.open(filename, mode='rt') as file:
        for record in SeqIO.parse(file, 'fasta'):
            records.append(record)
    return records


def main():


if __name__ == '__main__':
    main()
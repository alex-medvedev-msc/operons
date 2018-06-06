import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas
from os import path, makedirs, listdir

"""
feature	class	        assembly	    assembly_unit seq_type  chromosome	genomic_accession start	end	strand	product_accession	non-redundant_refseq	related_accession	name	symbol	GeneID	locus_tag	feature_interval_length	product_length	attributes
gene	protein_coding	GCF_000006665.1	Primary       Assembly	chromosome	NC_002655.2	      190	273	+		Z_RS00005	        84		old_locus_tag=Z0001
"""


def save_noncoding_regions(directory, filename, regions):
    regions_name = '{0}_regions.fasta'.format(filename)
    with open(path.join(directory, regions_name), 'w') as file:
        SeqIO.write(regions, file, 'fasta')


def extract_noncoding_regions(records, table: pandas.DataFrame):
    chromosome = records[0]
    if 'plasmid' in table["seq_type"].values:
        plasmid = records[1]
        if len(records) > 2:
            raise Exception('records have len({0}) instead of 2'.format(len(records)))

    regions = []
    genes = table[(table["# feature"] == 'CDS') & (table["class"] == 'with_protein')]
    coding_strand = genes[genes["strand"] == '+']
    template_strand = genes[genes["strand"] == '-']
    unique_id = 1

    for i in range(1, len(coding_strand)):

        end = coding_strand.iat[i-1, 8]
        start = coding_strand.iat[i, 7]
        seq_type = coding_strand.iat[i, 4]
        if seq_type == "plasmid":
            region = SeqRecord(plasmid.seq[end:start], id='{0}.{1}'.format(plasmid.id, unique_id))
        else:
            region = SeqRecord(chromosome.seq[end:start], id='{0}.{1}'.format(chromosome.id, unique_id))
        unique_id += 1
        regions.append(region)

    for i in range(1, len(template_strand)):

        end = template_strand.iat[i-1, 8]
        start = template_strand.iat[i, 7]
        seq_type = template_strand.iat[i, 4]

        if seq_type == "plasmid":
            region = SeqRecord(plasmid.seq[end:start].reverse_complement(), id='{0}.{1}'.format(plasmid.id, unique_id))
        else:
            region = SeqRecord(chromosome.seq[end:start].reverse_complement(), id='{0}.{1}'.format(chromosome.id, unique_id))
        unique_id += 1
        regions.append(region)

    return regions


def save_proteins(directory, filename, proteins):
    protein_name = '{0}_proteins.fasta'.format(filename)
    with open(path.join(directory, protein_name), 'w') as file:
        SeqIO.write(proteins, file, 'fasta')


def translate_record(source: SeqRecord, start, end, strand, id, name):
    if abs(end - start) > 30000:
        print('record {0} is too long: {1}'.format(source.id, abs(end - start)))
        return

    if strand == "+":
        protein = source.seq[start-1:end+1].translate()
    else:
        s3 = source.seq[start:end].reverse_complement()
        protein = s3.translate()

    if not isinstance(name, str):
        name = 'unknown protein'
    return SeqRecord(protein, id=id, name=name, description=name)


def translate(records, table: pandas.DataFrame):
    chromosome = records[0]
    if 'plasmid' in table["seq_type"].values:
        plasmid = records[1]
        if len(records) > 2:
            raise Exception('records have len({0}) instead of 2'.format(len(records)))

    proteins = []
    genes = table[(table["# feature"] == 'CDS') & (table["class"] == 'with_protein')]
    plasmid_index, chromosome_index = 0, 0
    for _, gene in genes.iterrows():
        id = gene["non-redundant_refseq"]
        if gene["seq_type"] == 'plasmid':
            if not isinstance(gene["non-redundant_refseq"], str):
                plasmid_index += 1
                id = '{0}_{1}'.format(plasmid.id, plasmid_index)

            protein = translate_record(plasmid, gene["start"], gene["end"], gene["strand"], id, gene["name"])
        else:
            if not isinstance(gene["non-redundant_refseq"], str):
                chromosome_index += 1
                id = '{0}_{1}'.format(chromosome.id, chromosome_index)

            protein = translate_record(chromosome, gene["start"], gene["end"], gene["strand"], id, gene["name"])

        if protein is not None:
            proteins.append(protein)

    return proteins


def load_feature_table(filename):
    ft_name = path.join('data', '{0}_feature_table.txt.gz'.format(filename))
    with gzip.open(ft_name, mode='rt') as file:
        table = pandas.read_csv(file, sep='\t')
        return table


def ungzip(filename):
    genomic_name = path.join('data', '{0}_genomic.fna.gz'.format(filename))
    records = []
    with gzip.open(genomic_name, mode='rt') as file:
        for record in SeqIO.parse(file, 'fasta'):
            records.append(record)
    return records


def get_genome_name(file):
    return file.replace('_genomic.fna.gz', '')


def count_full_genomes():
    full_count = 0
    others_count = 0
    for file in listdir('data'):
        if file.endswith('genomic.fna.gz'):
            genome_name = get_genome_name(file)
            records = ungzip(genome_name)
            if len(records) > 2:
                others_count = others_count + 1
            else:
                full_count = full_count + 1
    print('count of full genomes is {0}, count of others is {1}'.format(full_count, others_count))


def transform(genome_name, dest_dir):
    records = ungzip(genome_name)
    if len(records) > 2:
        print('genome {0} contains {1} scaffolds'.format(genome_name, len(records)))
        return
    feature_table = load_feature_table(genome_name)
    proteins = translate(records, feature_table)
    regions = extract_noncoding_regions(records, feature_table)
    save_proteins(dest_dir, genome_name, proteins)
    print('processed {1} proteins for genome {0}'.format(genome_name, len(proteins)))
    save_noncoding_regions(dest_dir, genome_name, regions)
    print('processed {1} non-coding regions for genome {0}'.format(genome_name, len(regions)))


def main():
    dest_dir = 'transformed'
    makedirs(dest_dir, exist_ok=True)
    for file in listdir('data'):
        if file.endswith('genomic.fna.gz'):
            genome_name = get_genome_name(file)
            transform(genome_name, dest_dir)

if __name__ == '__main__':
    main()

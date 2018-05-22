from ftplib import FTP
from time import sleep

home_dir = '/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/'


def download_file(ftp, directory):

    ftp.cwd(directory)

    feature_table = '{0}_feature_table.txt.gz'.format(directory)
    genome = '{0}_genomic.fna.gz'.format(directory)

    with open('data/{0}'.format(feature_table), 'wb') as feature_file:
        ftp.retrbinary('RETR {0}'.format(feature_table), feature_file.write)

    with open('data/{0}'.format(genome), 'wb') as genome_file:
        ftp.retrbinary('RETR {0}'.format(genome), genome_file.write)

    print("genome {0} downloaded".format(directory))

    ftp.cwd(home_dir)
    sleep(2)


def main():
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd(home_dir)
    dirs = ftp.nlst()
    with open('current_genome.txt', 'w') as cg:
        for d in dirs:
            download_file(ftp, d)
            cg.write(d)

    ftp.quit()


if __name__ == '__main__':
    main()

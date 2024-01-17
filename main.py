from gooey import Gooey, GooeyParser
import subprocess
import os
import sys
import threading

@Gooey(dump_build_config=True, program_name="16S MAT",default_size=(1000, 800),timing_options={'show_time_remaining': True, 'hide_time_remaining_on_complete': False})
def main():
    desc = "   16S MAT is a Windows-based application developed for 16S rRNA gene sequencing analysis."

    file_help_msg = "Select current fastq.gz file path"

    my_cool_parser = GooeyParser(description=desc)

    my_cool_parser.add_argument(
        "Input files", help=file_help_msg, widget="FileSaver")

    my_cool_parser.add_argument(
        "Reference database file 1", help="Select the reference sequence file for sequence alignment.\n(e.g. silva.bacteria.fasta)", widget="FileSaver")

    my_cool_parser.add_argument(
        "Reference database file 2", help="Select the reference sequence file for the classification algorithm.\n(e.g. trainset9_032012.pds.fasta)", widget="FileSaver")

    my_cool_parser.add_argument(
        "Reference database file 3",help="Select the reference taxonomy file.\n(e.g. trainset9_032012.pds.tax)", widget="FileSaver")

    my_cool_parser.add_argument('-d', '--criteria', default=90,help="The threshold percentage for screening sequences(default: 90)",dest='Criteria',
                                type=int)

    my_cool_parser.add_argument('-l', '--maxhomop',default=8,help="The maximum number of consecutive homologous base(default: 8)",dest='Maxhomop',
                                type=int)

    my_cool_parser.add_argument('-k', '--label', default=0.03,help="Minimum relative abundance threshold for OTU presence(default: 0.03)", dest='Label',
                                type=float)

    my_cool_parser.add_argument('-m', '--cutoff', default=0.03,help="Threshold for sequence similarity(default: 0.03)",dest='Cutoff',
                                type=float)

    my_cool_parser.add_argument(
        "-e", "--error", action="store_true", help="Stop process on error (default: No)")

    args = my_cool_parser.parse_args()
    input_files = args.__dict__["Input files"]
    reference_file1 = args.__dict__["Reference database file 1"]
    reference_file2 = args.__dict__["Reference database file 2"]
    reference_file3 = args.__dict__["Reference database file 3"]
    input_dir = os.path.dirname(input_files)
    output_dir = os.path.join(input_dir, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    mothur_cmd = f'mothur "#make.file(inputdir={os.path.dirname(input_files)}, outputdir={output_dir}, prefix=stability, type=gz) ; ' \
                 f'make.contigs(file={os.path.join(output_dir,"stability.files")}, processors=8) ; ' \
                 f'summary.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.fasta")}) ; ' \
                 f'screen.seqs(fasta={os.path.join(output_dir,"stability.trim.contigs.fasta")}, group={os.path.join(output_dir,"stability.contigs.groups")}, maxambig=0, optimize=maxlength, criteria={args.Criteria}); ' \
                 f'unique.seqs(fasta={os.path.join(output_dir,"stability.trim.contigs.good.fasta")}); ' \
                 f'count.seqs(name={os.path.join(output_dir,"stability.trim.contigs.good.names")}, group={os.path.join(output_dir,"stability.contigs.good.groups")}); ' \
                 f'align.seqs(fasta={os.path.join(output_dir,"stability.trim.contigs.good.unique.fasta")}, reference={reference_file1}, processors=8);' \
                 f'summary.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.align")}, count={os.path.join(output_dir, "stability.trim.contigs.good.count_table")}) ; ' \
                 f'screen.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.align")}, count={os.path.join(output_dir, "stability.trim.contigs.good.count_table")},summary={os.path.join(output_dir, "stability.trim.contigs.good.unique.summary")},optimize=start-end-minlength, criteria={args.Criteria}, maxhomop={args.Maxhomop}, processors=8); ' \
                 f'filter.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.align")}, vertical=T, trump=.); ' \
                 f'unique.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.fasta")}, count={os.path.join(output_dir, "stability.trim.contigs.good.good.count_table")}); ' \
                 f'pre.cluster(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.fasta")},count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.count_table")}, diffs=2); ' \
                 f'chimera.vsearch(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta")},count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table")}, dereplicate=t) ; ' \
                 f'remove.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta")}, accnos={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos")}, name={os.path.join(output_dir,"stability.trim.contigs.good.names")}); ' \
                 f'remove.lineage(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta")}, count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}, taxonomy={reference_file3}, taxon=unknown); ' \
                 f'classify.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta")}, count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}, reference={reference_file2}, taxonomy={reference_file3}, cutoff=80); ' \
                 f'remove.lineage(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta")}, count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}, taxonomy={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy")}, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota); ' \
                 f'summary.tax(taxonomy={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy")},count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}); ' \
                 f'dist.seqs(fasta={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta")}, cutoff={args.Cutoff}); ' \
                 f'cluster(column={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist")}, count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}); ' \
                 f'make.shared(list={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list")}, count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}, label={args.Label}); ' \
                 f'classify.otu(list={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list")}, count={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")}, taxonomy={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy")}, label={args.Label});' \
                 f'make.biom(shared={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared")},constaxonomy={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")});' \
                 f'rarefaction.single(shared={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared")}); ' \
                 f'summary.single(shared={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared")}, calc=shannon-invsimpson-chao); ' \
                 f'dist.shared(shared={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared")}, calc=thetayc-braycurtis, subsample=t); ' \
                 f'pcoa(phylip={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.dist")}); ' \
                 f'nmds(phylip={os.path.join(output_dir, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.dist")}, mindim = 2, maxdim = 5) ; ' \

    process = subprocess.Popen(mothur_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    verbosity = my_cool_parser.add_mutually_exclusive_group()

    with open('output.log', 'wb') as output_file:
        for line in process.stdout:
            output_file.write(line)
            try:
                print(line.decode(sys.getfilesystemencoding()), end='')
            except UnicodeDecodeError:
                pass

    process.wait()

    stdout, stderr = process.communicate()

    try:
        output = stdout.decode('utf-8')
    except UnicodeDecodeError:
        output = stdout.decode('gbk')

    keep_files = [
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.braycurtis.0.03.lt.ave.dist',
        'stability.files',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.groups.rarefaction',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.groups.summary',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.dist',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.nmds.axes',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.nmds.stress',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary',
        'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.dist']

    for file in os.listdir(output_dir):
        if file in keep_files or file.endswith('.svg'):
            continue
        file_path = os.path.join(output_dir, file)
        if os.path.isfile(file_path):
            os.remove(file_path)

    for file_name in keep_files:
        if file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.braycurtis.0.03.lt.ave.dist':
            new_file_name = 'Braycurtis_dist.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.groups.rarefaction':
            new_file_name = 'Rarefaction.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.groups.summary':
            new_file_name = 'Diversity.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared':
            new_file_name = 'Otu.table'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.dist':
            new_file_name = 'Thetayc_dist.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.nmds.axes':
            new_file_name = 'Nmds_axes.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.nmds.stress':
            new_file_name = 'Nmds_stress.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes':
            new_file_name = 'PCoA_axes.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary':
            new_file_name = 'Tax_summary.txt'
        elif file_name == 'stability.files':
            new_file_name = 'Sample_files.txt'
        elif file_name == 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.dist':
            new_file_name = 'Phylip_dist.txt'
        old_file_path = os.path.join(output_dir, file_name)
        new_file_path = os.path.join(output_dir, new_file_name)
        if os.path.exists(old_file_path):
            os.rename(old_file_path, new_file_path)

def run_main():
    t = threading.Thread(target=main)
    t.start()

if __name__ == '__main__':
    main()
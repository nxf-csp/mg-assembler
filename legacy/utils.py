import os
import subprocess
import shutil
import logging

logger_utils = logging.getLogger("default")

def convert_sam_to_bam(sam_file, bam_file, threads, exclude_unmapped=True):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    samtools = "samtools"
    
    flag = "-F 4" if exclude_unmapped else ""
    cmd = f"samtools view -@ {threads} -Sb {flag} -o {bam_file} {sam_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def sort_bam(unsorted_bam_file, sorted_bam_file, threads, by_name=False):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    samtools = "samtools"

    n_flag = "-n" if by_name else ""
    cmd = f"samtools sort -@ {threads} {n_flag} -o {sorted_bam_file} {unsorted_bam_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def index_bam(bam_file, threads):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    samtools = "samtools"
    
    cmd = f"samtools index -@ {threads} {bam_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def get_alignment_flagstat(aln_file, report_file, report_format, threads):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    samtools = "samtools"
    
    cmd = f"samtools flagstat -@ {threads} -O {report_format} {aln_file} > {report_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def run_seqkit_stats(fasta_list=[], fasta_dir="", ext="fasta", output_file="./bins_seq_stats.tsv"):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    seqkit = "seqkit"
    
    if fasta_list:
        input_files = " ".join(fasta_list)
    elif fasta_dir:
        input_files = os.path.join(fasta_dir, f"*.{ext}")
    else:
        logger_utils.error("Neither fasta list nor dirrectory with fasta files is provided!")
        return 1

    cmd = f"{seqkit} stats -T -a -N 90 {input_files} > {output_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def run_seqkit_replace(input_file, output_file, pattern="^", replacement=""):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    seqkit="seqkit"
    
    cmd = f"{seqkit} replace -p '{pattern}' -r '{replacement}' {input_file} >> {output_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def run_seqkit_grep(input_file, patterns_file, output_file):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    seqkit="seqkit"

    cmd = f"{seqkit} grep -f {patterns_file} {input_file} >> {output_file}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def run_bbmap_repair(r1_input, r2_input, r1_output, r2_output):
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    repair="repair.sh"

    cmd = f"""{repair} in={r1_input} in2={r2_input} out={r1_output} out2={r2_output} repair=t"""
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def remove_tmp(tmp_path):
    if os.path.isdir(tmp_path):
        try:
            shutil.rmtree(tmp_path)
            logger_utils.info(f"Directory '{tmp_path}' and its contents removed successfully.")
        except Exception as error:
            logger_utils.error(f"Unable to remove directory '{tmp_path}'.", exc_info=True)
    else:
        try:
            os.remove(tmp_path)
            logger_utils.info(f"File '{tmp_path}' removed successfully.")
        except Exception as error:
            logger_utils.error(f"Unable to remove file '{tmp_path}'.", exc_info=True)

#############################

def run_bowtie2_build(fasta_file, output_dir, threads):
    '''
    Запускает функцию bowtie2-build пакета Bowtie2. Создает индекс последовательности.

    Аргументы:
        fasta_file (str)  путь к FASTA файлу с нуклеотидной последовательсностью.
        output_dir (str)  путь к дирректории, куда будет сохранен индекс.
        threads (int)     число потоков.
    '''
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    bowtie2_build = "bowtie2-build"
    
    index_path = os.path.join(output_dir, "index")
    cmd = f"{bowtie2_build} --threads {threads} --quiet {fasta_file} {index_path}"
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1

def run_bowtie2(reads_dir, index_path, output_file, threads, use_unpaired=True):
    '''
    Запускает функцию bowtie2 пакета Bowtie2. Производит выравнивание ридов.

    Аргументы:
        reads_dir (str)      путь к директории с .fastq файлами предобработанных ридов.
        index_path (str)     путь к индексу в формате /path/to/dir/index.
        output_file (str)    путь к файлу выравнивания в формате .sam.
        threads (int)        число потоков.
        use_unpaired (bool)  использовать риды, утратившие пару в процессе тримминга? [True]
    '''
    env = "utils"
    conda_run_cmd = f"conda run -n {env}"
    bowtie2 = "bowtie2"

    unpaired_reads = ','.join([os.path.join(reads_dir, "R1_unpaired.fastq.gz"),
                               os.path.join(reads_dir, "R2_unpaired.fastq.gz")])
    u_flag = f"-U {unpaired_reads}" if use_unpaired else ""

    cmd = f"""{bowtie2} --threads {threads} -q --fr -x {index_path} 
              -1 {os.path.join(reads_dir, "R1_paired.fastq.gz")} 
              -2 {os.path.join(reads_dir, "R2_paired.fastq.gz")} {u_flag} 
              -S {output_file}"""
    cmd = " ".join(cmd.split())
    logger_utils.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_utils.error("An error occurred!", exc_info=True)
        return 1
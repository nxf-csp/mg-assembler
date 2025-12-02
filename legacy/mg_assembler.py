import os
import subprocess
import logging
import shutil
import json
from Bio import SeqIO
import pandas as pd

import utils

logger_asmb = logging.getLogger("default")

def merge_lanes(fastq_paths, megred_fastq):
    '''
    Объединяет несколько .fastq файлов.
    Аргументы:
        fastq_paths (list): список путей к .fastq файлам.
        megred_fastq (str): финальный .fastq файл.
    '''
    cmd = f"cat {' '.join(fastq_paths)} > {megred_fastq}"
    logger_asmb.debug(f"CMD: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_asmb.error("An error occurred!", exc_info=True)
        return 1

def trim_reads(r1_fastq, r2_fastq, output_dir, threads,
               report_title="fastp_report", report_dir="",
               min_q=20, min_len=70, deduplication=True):
    '''
    Запускает fastp.
    В output_dir соханяются R1_paired.fastq.gz и R2_paired.fastq.gz с парными ридами,
    R1_unpaired.fastq.gz и R2_unpaired.fastq.gz с потерявшими пару ридами и
    discarded.fastq.gz с отсеяннами ридами.
    В report_dir сохраняются репорты fastp_report.json и fastp_report.html.
    По умолчанию выполняется автодетекция адаптеров для парных ридов и 
    удаление гомополимерных последовательностей.
    
    Аргументы:
        r1_fastq (str)        путь к .fastq файлу с R1 ридами.
        r2_fastq (str)        путь к .fastq файлу с R2 ридами.
        output_dir (str)      дирректория для сохранения результатов тримминга.
        threads (int)         число потоков.
        report_title (str)    заголовок репортов.
        report_dir (str)      дирректория для сохранения репортов [""]. 
                              Если report_dir="", то репорты сохраняются в output_dir.
        min_q (int)           Минимальное качество основания [20].
        min_len (int)         Минимальная допустимая длина рида после тримминга [70].
        deduplication (bool)  Дедуплицировать риды? [True] Дедупликация осуществляется со
                              стандартными утсановками fastp.
    '''
    env = "assembly"
    conda_run_cmd = f"conda run -n {env}"
    fastp = "fastp"
    
    dedup = "--dedup" if deduplication else ""
    report_dir = report_dir if report_dir else output_dir
    cmd = f"""{fastp} --thread {threads} 
              -q {min_q} -l {min_len} {dedup} --trim_poly_x --detect_adapter_for_pe 
              --unpaired1 {os.path.join(output_dir, "R1_unpaired.fastq.gz")} 
              --unpaired2 {os.path.join(output_dir, "R2_unpaired.fastq.gz")} 
              --failed_out {os.path.join(output_dir, "discarded.fastq.gz")} 
              --json {os.path.join(report_dir, "fastp_report.json")} 
              --html {os.path.join(report_dir, "fastp_report.html")} 
              --report_title '{report_title}' 
              --in1 {r1_fastq} 
              --in2 {r2_fastq} 
              --out1 {os.path.join(output_dir, "R1_paired.fastq.gz")} 
              --out2 {os.path.join(output_dir, "R2_paired.fastq.gz")}"""
    cmd = " ".join(cmd.split())
    logger_asmb.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_asmb.error("An error occurred!", exc_info=True)
        return 1

def run_metaspades(trimmed_reads_dir, output_dir, threads, memory_limit):
    '''
    Запускает MetaSPAdes.
    В output_dir соханяются все результаты сборки метагенома.
    Сборка осуществляется с использованием как парных, так и потерявших пару ридов и
    стандартным набором k-меров [21, 33, 55].
    
    Аргументы:
        trimmed_reads_dir (str)  путь к директории с предобработанными .fastq файлами.
        output_dir (str)         дирректория для сохранения результатов сборки.
        threads (int)            число потоков.
        memory_limit (int)       выделяемый объем памяти в Gb.
    '''
    env = "assembly"
    conda_run_cmd = f"conda run -n {env}"
    metaspades = "metaspades.py"

    cmd = f"""{metaspades} -t {threads} -m {memory_limit} 
              --pe1-1 {os.path.join(trimmed_reads_dir, 'R1_paired.fastq.gz')} 
              --pe1-2 {os.path.join(trimmed_reads_dir, 'R2_paired.fastq.gz')} 
              --pe1-s {os.path.join(trimmed_reads_dir, 'R1_unpaired.fastq.gz')} 
              --pe1-s {os.path.join(trimmed_reads_dir, 'R2_unpaired.fastq.gz')} 
              -o {output_dir}"""
    cmd = " ".join(cmd.split())
    logger_asmb.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_asmb.error("An error occurred!", exc_info=True)
        return 1

def process_contigs(unfiltered_contigs, filtered_contigs, 
                    min_len=1000, contig_basename="contig", 
                    no_description=False, add_to_description=""):
    '''
    Осуществляет постобработку контигов: фильтрует по длине, переименовывает и 
    добавляет описание. Описание формируется из длины контига, покрытия k-мерами и
    дополнительной информации, переданной в add_to_description (опционально).

    Использует SeqIO из Bio
    
    Аргументы:
        unfiltered_contigs (str)  путь к оригинальному .fasta файлу, генерируемому MetaSPAdes.
        filtered_contigs (str)    путь к .fasta файлу с обработанными контигами.
        min_len (int)             минимальная длина контига [1000].
        contig_basename (str)     основа имени контига ["contig"].
        no_description (bool)     не добавлять описание к записям итогового .fasta файла? [False]
        add_to_description (str)  доп. информация в описание [""]
    '''
    try:
        with open(unfiltered_contigs) as i_file:
            with open(filtered_contigs, 'a') as o_file:
                for contig in SeqIO.parse(i_file, "fasta"):
                    if len(contig.seq) < min_len:
                        continue
                    contig_name = contig.id.split("_")
                    contig.id = f"{contig_basename}_{contig_name[1]}" 
                    if no_description:
                        contig.description = ""
                    else:
                        contig.description = f"length={len(contig.seq)};spades_cov={contig_name[5]};{add_to_description}"
                    SeqIO.write(contig, o_file, "fasta")
        return 0
    except Exception as error:
        logger_asmb.error("An error occurred!", exc_info=True)
        return 1

def assemble_mg(sample, threads, memory_limit, 
                min_q=20, read_min_len=70, deduplication=True,
                contig_min_len=1000, contig_basename="contig", 
                clean_up=True):
    '''
    Реализует полный процесс сборки: 
    0)  Проверяет существование .fastq файлов с ридами.
    1)  fastq файлы копируются в директорию сборки образца. Если для R1 и R2 
        получено более одного .fastq файла, они объединяются в общий файл.
    2)  С использованием .fastp осуществляется тримминг и фильрация ридов. 
        Обработанные риды сохраняются в директории "reads".
    3)  Прошедшие тримминг парные и потерявшие пару риды передаются MetaSPAdes для сборки метагенома. 
    4)  Полученные контиги фильтруются, переименовывыются и сохраняются в директории "assembly".
    5)  Другие результаты сборки, представляющие интерес: граф сборки, пути контигов и каффолдов, 
        сохраняются в директории "assembly".
    6)  Удаляет временные файлы, копии необработанных .fastq файлов с ридами и другие результаты сборки.
    В процессе работы собирается статистика тримминга и полученнйо сборки.

    В результате работы в директории сборки образца остаются:
    sample_dir
    |----- assembly.log
    |----- assembly
    |      |----- contigs.fasta
    |      |----- assembly_graph.fastg
    |      |----- contigs.paths
    |      |----- scaffolds.paths
    |----- reads
    |      |----- R1_paired.fastq.gz
    |      |----- R2_paired.fastq.gz
    |      |----- R1_unpaired.fastq.gz
    |      |----- R2_unpaired.fastq.gz
    |      |----- discarded.fastq.gz
    |----- reports
           |----- spades.log
           |----- fastp_report.html
           |----- fastp_report.json
    
    Использует:
        - pandas
        - shutil
        - json
        - sample_cl.py
        - utils.py
        
    Аргументы:
        sample                    экземпляр класса AssembledSample, созданный для конкретного образца,
                                  определяющий структуру содержания директории сборки образца и 
                                  хранящий информацию о статусе сборки, файлах с важными результатами,
                                  а также собираемые метрики.                                  
        threads (int)             число потоков.
        memory_limit (int)        выделяемый объем памяти в Gb.
        min_q (int)               минимальное качество основания [20].
        min_len (int)             минимальная допустимая длина рида после тримминга [70].
        deduplication (bool)      дедуплицировать риды? [True] Дедупликация осуществляется со
                                  стандартными утсановками fastp.
        min_len (int)             минимальная длина контига, п.о. [1000].
        contig_basename (str)     основа имени контига ["contig"].
        clean_up (bool)           удалить все промежуточные файлы? [True] 

    Возвращает:
        - статус сборки образца ["Assembled" | "Assembled, warnings" | "Assembly failed"]
        - словарь с информацией о статусах подпроцессов, файлах с результатами и метриками.
    '''
    logger_asmb.debug(f"Working in '{os.getcwd()}'.")

    logger_asmb.info("Verifying inputs...")
    try:
        sample.check_inputs()
    except FileNotFoundError as error:
        logger_asmb.error("Inputs verification FAILED!", exc_info=True)
        logger_asmb.critical("Terminating sample processing.")
        sample.set_process_status("inputsCheck", "FAILED")
        sample.status = "Assembly failed"
        return sample
    sample.set_process_status("inputsCheck", "COMPLETED")
    logger_asmb.info("Inputs verification COMPLETED.")
    
    logger_asmb.info("Preprocessing .fastq files...")
    os.makedirs(sample.reads)
    os.makedirs(sample.reports)
    for r in sample.fastqs.keys():
        megred_fastq = os.path.join(sample.reads, f"{r}.fastq.gz")
        if len(sample.fastqs[r]) > 1:
            logger_asmb.info(f"Several {r} .fastq files recieved.")
            task_status = merge_lanes(fastq_paths=sample.fastqs[r], 
                                      megred_fastq=megred_fastq)
            if task_status:
                logger_asmb.critical(f"Merging of {r} .fastq files FAILED!")
                logger_asmb.critical("Terminating sample processing.")
                sample.set_process_status("readsProcessing", "FAILED")
                sample.status = "Assembly failed"
                return sample
            logger_asmb.info(f"Received {r} .fastq files are merged in '{megred_fastq}'.")
        else:
            try:
                shutil.copy(sample.fastqs[r][0], megred_fastq)
                logger_asmb.info(f"Received {r} .fastq file is copied to '{megred_fastq}'.")
            except Exception as error:
                logger_asmb.error("An error occurred!", exc_info=True)
                logger_asmb.critical(f"Copying {r} .fastq file FAILED!")
                logger_asmb.critical("Terminating sample processing.")
                sample.set_process_status("readsProcessing", "FAILED")
                sample.status = "Assembly failed"
                return sample
    logger_asmb.info("Preprocessing .fastq files COMPLETED.")

    logger_asmb.info("Starting reads trimming...")
    task_status = trim_reads(r1_fastq=os.path.join(sample.reads, "R1.fastq.gz"), 
                             r2_fastq=os.path.join(sample.reads, "R2.fastq.gz"), 
                             output_dir=sample.reads,
                             threads=threads,
                             report_dir=sample.reports,
                             report_title=f"{sample.pool}:{sample.id}",
                             min_q=min_q,
                             min_len=read_min_len, 
                             deduplication=deduplication)
    if task_status:
        logger_asmb.critical("Trimming of reads FAILED!")
        logger_asmb.critical("Terminating sample processing.")
        sample.set_process_status("readsProcessing", "FAILED")
        sample.status = "Assembly failed"
        return sample
    
    logger_asmb.info("Collecting trimming statistics...")
    try:
        with open(os.path.join(sample.reports, "fastp_report.json"), "r") as report:
            trimming_stats = json.load(report)["summary"]
        sample.memoize_metric("readsBeforeTrimming", trimming_stats['before_filtering']['total_reads'])
        sample.memoize_metric("readsAfterTrimming", trimming_stats['after_filtering']['total_reads'])
        r_prcnt = round(100 * trimming_stats['after_filtering']['total_reads'] / trimming_stats['before_filtering']['total_reads'], 2)
        sample.memoize_metric("prcntReadsPassed", r_prcnt)
    except Exception as error:
        logger_asmb.error("An error occurred!", exc_info=True)
    sample.set_process_status("readsProcessing", "COMPLETED")
    logger_asmb.info("Trimming of reads COMPLETED.")
    
    logger_asmb.info("Starting assembling...")
    os.makedirs(sample.swd)
    task_status = run_metaspades(trimmed_reads_dir=sample.reads, 
                                 output_dir=sample.swd, 
                                 threads=threads, 
                                 memory_limit=memory_limit)
    if task_status:
        logger_asmb.critical("Assembling FAILED!")
        logger_asmb.critical("Terminating sample processing.")
        sample.set_process_status("assembling", "FAILED")
        sample.status = "Assembly failed"
        return sample
    sample.set_process_status("assembling", "COMPLETED")
    logger_asmb.info("Assembling COMPLETED.")

    logger_asmb.info("Processing assembly results...")
    os.makedirs(sample.assembly)
    got_warnings = False
        
    logger_asmb.info("Filtering and renaming contigs...")
    task_status = process_contigs(unfiltered_contigs=os.path.join(sample.swd, "contigs.fasta"), 
                                  filtered_contigs=sample.outputs["assemblyFasta"], 
                                  min_len=contig_min_len, 
                                  contig_basename=contig_basename, 
                                  no_description=True)
    if task_status:
        logger_asmb.warning("Unable to finish contigs processing.")
        got_warnings = True

    logger_asmb.info("Collecting assembly metrics...")
    assembly_metrics_tsv = os.path.join(sample.swd, "assembly_metrics.tsv")
    task_status = utils.run_seqkit_stats(fasta_list=[sample.outputs["assemblyFasta"]], 
                                         output_file=assembly_metrics_tsv)
    if task_status:
        logger_asmb.warning("Unable to collect assembly metrics.")
        got_warnings = True
    else:
        try:
            assembly_metrics = pd.read_csv(assembly_metrics_tsv, sep="\t")
            sample.memoize_metric("contigsNum", int(assembly_metrics.at[0, "num_seqs"]))
            sample.memoize_metric("contigsSumLength", int(assembly_metrics.at[0, "sum_len"]))
            sample.memoize_metric("contigMinLength", int(assembly_metrics.at[0, "min_len"]))
            sample.memoize_metric("contigMaxLength", int(assembly_metrics.at[0, "max_len"]))
            sample.memoize_metric("NsNum", int(assembly_metrics.at[0, "sum_n"]))
        except Exception as error:
            logger_asmb.error("An error occurred!", exc_info=True)
    
    logger_asmb.info("Moving assembly graph file...")
    try:
        shutil.move(os.path.join(sample.swd, "assembly_graph.fastg"),
                    sample.outputs["assemblyGraph"])
    except Exception as error:
        logger_asmb.error("Unable to move assembly graph file.", exc_info=True)
        got_warnings = True
    logger_asmb.info("Moving contigs and scaffolds paths files...")
    try:
        shutil.move(os.path.join(sample.swd, "contigs.paths"),
                    sample.outputs["contigPaths"])
        shutil.move(os.path.join(sample.swd, "scaffolds.paths"),
                    sample.outputs["scaffoldPaths"])
    except Exception as error:
        logger_asmb.error("Unable to move contigs and scaffolds paths files.", exc_info=True)
        got_warnings = True
    logger_asmb.info("Saving MetaSPAdes log file...")
    try:
        shutil.move(os.path.join(sample.swd, "spades.log"),
                    os.path.join(sample.reports, "spades.log"))
        logger_asmb.info(f"MetaSPAdes log file has been moved to '{sample.reports}'")
    except Exception as error:
        logger_asmb.error("Unable to move metaSPAdes log file.", exc_info=True)

    if got_warnings:
        logger_asmb.warning("Processing assembly results FINISHED WITH WARNINGS!")
        logger_asmb.info(f"Temporary files of the assembling will be saved in '{sample.swd}'")
        sample.set_process_status("assemblyProcessing", "FAILED")
        sample.status = "Assembled (warnings)"
        return sample
    sample.set_process_status("assemblyProcessing", "COMPLETED")
    sample.status = "Assembled"

    if clean_up:
        logger_asmb.info("Removing temporary files...")
        tmps = [sample.swd, 
                os.path.join(sample.reads, "R1.fastq.gz"),
                os.path.join(sample.reads, "R2.fastq.gz")]
        for tmp in tmps:
            utils.remove_tmp(tmp)
        
    return sample
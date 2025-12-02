import os
import subprocess
import concurrent.futures
import threading
import logging
import shutil
from Bio import SeqIO
import pandas as pd

import utils

logger_ant = logging.getLogger("default")

##### Таксономическая идентификация и QC бинов #####

def identify_gtdbtk(genome_dir, work_dir, threads, 
                    tmp_dir="/tmp", ext="fasta"):
    '''
    Запускает функцию identify пайплайна GTDB-Tk.
    Аргументы:
        genome_dir (str)   путь к директории с последовательностями геномов.
        work_dir (str)     путь к дирректории, куда будут сохраняться результаты работы GTDB-Tk.
                           В ней будет создана директория "identify", куда будут записаны результаты
                           работы данной функции.
        tmp_dir (str)      дирректория для временных файлов.
        threads (int)      число потоков.
        ext (str)          расширение FASTA файлов ["fasta"]. 
        
    '''
    env = "gtdbtk"
    conda_run_cmd = f"conda run -n {env}"
    gtdbtk="gtdbtk"
    
    cmd = f"""{gtdbtk} identify 
               --cpus {threads} 
               --genome_dir {genome_dir} 
               --extension {ext} 
               --out_dir {os.path.join(work_dir, "identify")} 
               --tmpdir {tmp_dir}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def align_gtdbtk(work_dir, threads, tmp_dir="/tmp"):
    '''
    Запускает функцию align пайплайна GTDB-Tk.
    Аргументы:
        work_dir (str)     путь к дирректории, куда будут сохраняться результаты работы GTDB-Tk.
                           В ней будет создана директория "align", куда будут записаны результаты
                           работы данной функции.
        tmp_dir (str)      дирректория для временных файлов.
        threads (int)      число потоков.
        
    '''
    env = "gtdbtk"
    conda_run_cmd = f"conda run -n {env}"
    gtdbtk="gtdbtk"
    
    cmd = f"""{gtdbtk} align 
               --cpus {threads} 
               --identify_dir {os.path.join(work_dir, "identify")} 
               --out_dir {os.path.join(work_dir, "align")} 
               --tmpdir {tmp_dir}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def classify_gtdbtk(work_dir, genome_dir, threads, tmp_dir="/tmp", ext="fasta"):
    '''
    Запускает функцию classify пайплайна GTDB-Tk.
    Аргументы:
        work_dir (str)     путь к дирректории, куда будут сохраняться результаты работы GTDB-Tk.
                           В ней будет создана директория "classify", куда будут записаны результаты
                           работы данной функции.
        genome_dir (str)   путь к директории с последовательностями геномов.
        tmp_dir (str)      дирректория для временных файлов.
        threads (int)      число потоков.
        ext (str)          расширение FASTA файлов ["fasta"]. 
        
    '''
    env = "gtdbtk"
    conda_run_cmd = f"conda run -n {env}"
    gtdbtk="gtdbtk"
    
    cmd = f"""{gtdbtk} classify 
               --cpus {threads} 
               --genome_dir {genome_dir} 
               --extension {ext} 
               --align_dir {os.path.join(work_dir, "align")} 
               --out_dir {os.path.join(work_dir, "classify")}  
               --tmpdir {tmp_dir} 
               --skip_ani_screen"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_gtdbtk(bins_dir, output_dir, threads, bin_extension="fasta"):
    '''
    Последовательно запускает все этапы пайплайна GTDB-Tk.
    Аргументы:
        bins_dir (str)       путь к директории с последовательностями геномов (бинами).
        output_dir (str)     путь к дирректории, куда будут сохраняться результаты работы GTDB-Tk.
        threads (int)        число потоков.
        bin_extension (str)  расширение FASTA файлов бинов ["fasta"]. 
        
    '''
    os.makedirs(output_dir, exist_ok=True)
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    
    logger_ant.info("Gene calling...")
    task_status = identify_gtdbtk(genome_dir=bins_dir, 
                                  work_dir=output_dir, 
                                  threads=threads,
                                  ext=bin_extension)
    if task_status:
        logger_ant.critical("Gene calling FAILED!")
        return 1
    logger_ant.info("Gene calling COMPLETED.")
    
    logger_ant.info("Alignment of marker genes sequences...")
    task_status = align_gtdbtk(work_dir=output_dir, 
                               threads=threads)
    if task_status:
        logger_ant.critical("Alignment of marker genes sequences FAILED!")
        return 1
    logger_ant.info("Alignment of marker genes sequences COMPLETED.")
        
    logger_ant.info("Classifying bins...")
    task_status = classify_gtdbtk(work_dir=output_dir, 
                                  genome_dir=bins_dir,
                                  threads=threads,
                                  ext=bin_extension)
    if task_status:
        logger_ant.critical("Bins classification FAILED!")
        return 1
    logger_ant.info("Bins classification COMPLETED.")
    return 0

def move_gtdbtk_gene_calls(gtdbtk_output_dir, dest_dir, ext="faa"):
    '''
    Перемещает результаты предсказния генов, полученные в процессе работы GTDB-Tk,
    в общую директорию. Написана специально под структуру результатов GTDB-Tk.

    Использует:
        - shutil
    
    Аргументы:
        gtdbtk_output_dir (str)  путь к директории с результатами работы GTDB-Tk.
        dest_dir (str)           путь к дирректории, куда будут перемещены файлы.
        ext (str)                расширение файлов: .faa, .fna, .gff ["faa"]. 
        
    '''
    bins_gc_dir = f"{gtdbtk_output_dir}/identify/identify/intermediate_results/marker_genes"
    for mg_bin in os.listdir(bins_gc_dir):
        try:
            shutil.copyfile(os.path.join(bins_gc_dir, mg_bin, f"{mg_bin}_protein.{ext}"), 
                            os.path.join(dest_dir, f"{mg_bin}.genes.{ext}"))
        except Exception:
            logger_ant.error(f"{mg_bin}: an error occurred!", exc_info=True)
    if len(os.listdir(bins_gc_dir)) == len(os.listdir(dest_dir)):
        logger_ant.info("All gene call files were copied successfully.")
    else:
        logger_ant.warning("Copying some gene call files FAILED!")

def run_checkm(bins_dir, output_dir, threads, ext="fasta"):
    '''
    Запускает CheckM для анализа качества геномов (бинов).
    Аргументы:
        bins_dir (str)     путь к директории с последовательностями геномов.
        output_dir (str)   путь к дирректории, куда будут сохраняться результаты работы CheckM.
        threads (int)      число потоков.
        ext (str)          расширение FASTA файлов ["fasta"]. Если ext="faa", QC осуществляется на 
                           основе переданных аминокислотных последовательностей CDS; 
                           CheckM не осуществляет предсказание генов.
    '''
    env = "checkm"
    conda_run_cmd = f"conda run -n {env}"
    checkm = "checkm"

    genes = "--genes" if ext == "faa" else ""
        
    cmd = f"""{checkm} lineage_wf 
              -t {threads} 
              {genes} -x {ext} 
              --tab_table -f {os.path.join(output_dir, "checkm_report.tsv")} 
              {bins_dir} {output_dir}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_checkm2(bins_dir, output_dir, threads, ext="fasta"):
    '''
    Запускает CheckM2 для анализа качества геномов (бинов).
    Аргументы:
        bins_dir (str)     путь к директории с последовательностями геномов.
        output_dir (str)   путь к дирректории, куда будут сохраняться результаты работы CheckM2.
        threads (int)      число потоков.
        ext (str)          расширение FASTA файлов ["fasta"]. Если ext="faa", QC осуществляется на 
                           основе переданных аминокислотных последовательностей CDS; 
                           CheckM2 не осуществляет предсказание генов.
    '''
    env = "checkm2"
    conda_run_cmd = f"conda run -n {env}"
    checkm2 = "checkm2"

    genes = "--genes" if ext == "faa" else ""
        
    cmd = f"""{checkm2} predict 
              -t {threads} 
              {genes} -x {ext} 
              -i {bins_dir} -o {output_dir}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def bins_qc(bins_dir, output_dir, threads, bin_extension="fasta"):
    '''
    Запускает CheckM и CheckM2 в параллельных потоках. Выделенные потоки CPU
    распределяются в соотношении 11:1 между CheckM и CheckM2 соотвественно.
    
    Использует:
        - concurrent.futures
        - threading
    
    Аргументы:
        bins_dir (str)     путь к директории с последовательностями геномов.
        output_dir (str)   путь к дирректории, куда будут сохраняться результаты работы CheckM и CheckM2.
        threads (int)      число потоков.
        ext (str)          расширение FASTA файлов ["fasta"].
    '''
    
    checkm2_threads = round(threads / 12)
    checkm_threads = threads - checkm2_threads - 1
    
    def run_checkm_task():
        logger_ant.info("Bins QC with CheckM...")
        task_status = run_checkm(bins_dir=bins_dir, 
                                 output_dir=os.path.join(output_dir, "checkm"),  
                                 threads=checkm_threads,
                                 ext=bin_extension)
        if task_status:
            logger_ant.warning("Bins QC with CheckM FAILED!")
            return
        else:
            logger_ant.info("Bins QC with CheckM COMPLETED.")
            return

    def run_checkm2_task():
        logger_ant.info("Bins QC with CheckM2...")
        task_status = run_checkm2(bins_dir=bins_dir, 
                                  output_dir=os.path.join(output_dir, "checkm2"),  
                                  threads=checkm2_threads,
                                  ext=bin_extension)
        if task_status:
            logger_ant.warning("Bins QC with CheckM2 FAILED!")
            return
        else:
            logger_ant.info("Bins QC with CheckM2 COMPLETED.")
            return

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
    futures = [executor.submit(run_checkm_task),
               executor.submit(run_checkm2_task)]
    executor.shutdown(wait=True)

def merge_tax_and_qc(gtdbtk_output_dir, checkm_output_dir, checkm2_output_dir, 
                     bins_stats_tsv, output_file, checkm_v=1):
    '''
    Объединяет репорты GTDB-Tk (gtdbtk.bac120.summary.tsv и gtdbtk.ar53.summary.tsv),
    CheckM (checkm_report.tsv), CheckM2 (quality_report.tsv) и метрики MAGs (bins_stats_tsv) 
    в общий репорт. 
    Собирает информацию о числе MAGs разного качества (HQ, MQ и LQ). 
    Возвращает в виде словаря следующей структуры:
    MAGs = {"HQ" : #, "MQ" : #, "LQ"}
    
    Использует:
        - pandas
    
    Аргументы:
        gtdbtk_output_dir (str)   путь к дирректории c результатами работы GTDB-Tk.
        checkm_output_dir (str)   путь к дирректории c результатами работы CheckM.
        checkm2_output_dir (str)  путь к дирректории c результатами работы CheckM2.
        bins_stats_tsv (str)      путь к .tsv файлу с метриками MAGs (utils.run_seqkit_stats()).
        output_file (str)         путь к .tsv файлу, в который будет записан обобщенный репорт.
        checkm_v (int)            версия CheckM, чьи метрики завершенности и контаминации использовать
                                  для присвоения категории качества MAGs [1]. CheckM - 1; CheckM2 - 2.
    '''
    required_files = {"gtdbtk.bac120" : os.path.join(gtdbtk_output_dir, "classify", "gtdbtk.bac120.summary.tsv"),
                      "checkm" : os.path.join(checkm_output_dir, "checkm_report.tsv"),
                      "checkm2" : os.path.join(checkm2_output_dir, "quality_report.tsv"),
                      "bins_stats" : bins_stats_tsv}
    for f in required_files.keys():
        if os.path.exists(required_files[f]):
           logger_ant.info(f"File '{required_files[f]}' found.")
        else:
            logger_ant.error(f"The required '{required_files[f]}' does not exist. Terminating.")
            return 1, None

    gtdbtk_report = pd.read_csv(required_files["gtdbtk.bac120"], sep="\t")
    gtdbtk_ar53 = os.path.join(gtdbtk_output_dir, "classify", "gtdbtk.ar53.summary.tsv")
    if os.path.exists(gtdbtk_ar53):
        gtdbtk_report = pd.concat([gtdbtk_report, pd.read_csv(gtdbtk_ar53, sep="\t")])
        
    lf = lambda r: "Unclassified" if "Unclassified" in r["classification"] else r["classification"]
    gtdbtk_report["classification"] = gtdbtk_report.apply(lf, axis=1)
    gtdbtk_report["organism"] = gtdbtk_report["classification"]
    gtdbtk_report["organism"] = gtdbtk_report["organism"].apply(lambda r: r.split(";")[-1].split("__")[-1])
    gtdbtk_report["MAG_q"] = pd.NA
    gtdbtk_report.loc[gtdbtk_report["organism"] == "Unclassified", "organism"] = pd.NA
    gtdbtk_report = gtdbtk_report[["user_genome", "MAG_q", "organism", "classification", "closest_genome_reference"]]
    gtdbtk_report = gtdbtk_report.rename(columns = {"user_genome" : "bin"})
    
    checkm_report = pd.read_csv(required_files["checkm"], sep="\t")
    checkm_report = checkm_report[["Bin Id", "Completeness", "Contamination", "Strain heterogeneity"]]
    checkm_report = checkm_report.rename(columns={"Bin Id" : "bin", 
                                                  "Completeness" : "checkm_completeness", 
                                                  "Contamination" : "checkm_contamination", 
                                                  "Strain heterogeneity" : "checkm_strain_heterogeneity"})
    checkm_report["bin"] = checkm_report["bin"].apply(lambda r: r.split(".")[0])
    checkm_report.loc[len(checkm_report)] = ["bin_000", 0, 0, 0]
    
    checkm2_report = pd.read_csv(required_files["checkm2"], sep="\t")
    checkm2_report = checkm2_report[["Name", "Completeness", "Contamination"]]
    checkm2_report = checkm2_report.rename(columns={"Name" : "bin", 
                                                    "Completeness" : "checkm2_completeness", 
                                                    "Contamination" : "checkm2_contamination"})
    checkm2_report["bin"] = checkm2_report["bin"].apply(lambda r: r.split(".")[0])
    checkm2_report.loc[len(checkm2_report)] = ["bin_000", 0, 0]
    
    bins_seq_stats = pd.read_csv(required_files["bins_stats"], sep="\t")
    bins_seq_stats = bins_seq_stats[["file", "GC(%)", "num_seqs", "sum_len", "min_len", "max_len", 
                                     "N50", "N50_num", "N90", "sum_gap", "sum_n"]]
    bins_seq_stats = bins_seq_stats.rename(columns={"file" : "bin",
                                                    "GC(%)" : "GC",
                                                    "N50_num" : "L50"})
    bins_seq_stats["bin"] = bins_seq_stats["bin"].apply(lambda r: r.split("/")[2].split(".")[0])
    
    bins_stats = pd.merge(gtdbtk_report, checkm_report, how="right", on="bin")
    bins_stats = bins_stats.merge(checkm2_report, how="left", on="bin")
    bins_stats = bins_stats.merge(bins_seq_stats, how="left", on="bin")
    bins_stats[bins_stats.isna()] = pd.NA

    def estimate_mags_q(r):
        if checkm_v == 1:
            commpleteness = r["checkm_completeness"]
            contamination = r["checkm_contamination"]
        elif checkm_v == 2:
            commpleteness = r["checkm2_completeness"]
            contamination = r["checkm2_contamination"]

        if commpleteness >= 90 and contamination < 5:
            return "HQ"
        elif commpleteness >= 50 and contamination < 10:
            return "MQ"
        else:
            return "LQ"

    bins_stats["MAG_q"] = bins_stats.apply(estimate_mags_q, axis=1)
    bins_stats["MAG_q"] = pd.Categorical(bins_stats["MAG_q"], categories=["LQ", "MQ", "HQ"], ordered=True)
    if checkm_v == 1:
        use_completeness = "checkm_completeness"
    elif checkm_v == 2:
        use_completeness = "checkm2_completeness"
    bins_stats = bins_stats.sort_values(by=["MAG_q", use_completeness], ascending=False)
    
    MAGs = bins_stats["MAG_q"].value_counts().to_dict()
    MAGs = {k : MAGs[k] for k in ["HQ", "MQ", "LQ"]}
    
    try:
        bins_stats.to_csv(output_file, sep='\t', index=False)
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1, None
    return 0, MAGs

def bins_tax_ident_and_qc(bins_dir, work_dir, bins_metrics_tsv,
                          threads, bin_extension="fasta", 
                          use_aa_seqs_for_qc=True, checkm_v=1):
    '''
    Реализует этап таксономической идентификации и оценки качества полученных бинов,
    включающий следующие подэтапы:
    0)  Проверяет существование FASTA файлов с нужным расширением в указанной директории.
    1)  Собирает метрики последовательностей бинов с сипользованием функции run_seqkit_stats()
        модуля utils.
    2)  Осуществляет таксономическую идентификацию бинов и предсказание CDS в их последовательностях
        посредством пайплайна GTDB-Tk. Бин bin_000, включающий в себя все контиги, 
        которые не вошли в другие бины, маскируется и не анализируется GTDB-Tk.
    3)  Осуществляет QC бинов на основе предсказаных CDS посредством CheckM и CheckM2. 
        Бин bin_000, включающий в себя все контиги, которые не вошли в другие бины, 
        маскируется и не анализируется CheckM и CheckM2. 
    4)  Объединяет результаты работы GTDB-Tk, CheckM и CheckM2, а также метрики MAGs в общий репорт;
        демаскирует bin_000.
    В процессе работы собирает информацию о числе MAGs разного качества (HQ, MQ и LQ).

    Использует:
        - utils.py

    Аргументы:
        bins_dir (str)             путь к директории с последовательностями геномов.
        work_dir (str)             путь к дирректории, куда будут сохраняться результаты.
        bins_metrics_tsv (str)     путь к .tsv файлу, в который будет записана собранная информация о бинах.
        threads (int)              число потоков.
        bin_extension (str)        расширение FASTA файлов ["fasta"].
        use_aa_seqs_for_qc (bool)  Использовать аминокислотные последовательности CDS для QC бинов? Если нет,
                                   то QC бинов будет осуществляться на полных нуклеотидных последовательностях.
        checkm_v (int)             версия CheckM, чьи метрики завершенности и контаминации использовать
                                   для присвоения категории качества MAGs [1]. CheckM - 1; CheckM2 - 2.
    '''
    logger_ant.info("Checking bins...")
    if os.path.exists(bins_dir):
        n_bins = len([f for f in os.listdir(bins_dir) if f.endswith(bin_extension)])
        if n_bins:
            logger_ant.info(f"Directory with {n_bins} .{bin_extension} files found.")
        else:
            logger_ant.critical(f"There are no .{bin_extension} files in the specified directory!")
            logger_ant.critical("Terminating bins taxonomy identification and QC.")
            return 1, None
    else:
        logger_ant.critical(f"Directory with bins .{bin_extension} files not found!")
        logger_ant.critical("Terminating bins taxonomy identification and QC.")
        return 1, None

    logger_ant.info("Collecting sequence metrics for each bin...")
    bins_seq_stats = os.path.join(work_dir, "bins_seq_stats.tsv")
    task_status = utils.run_seqkit_stats(fasta_dir=bins_dir, 
                                         ext=bin_extension, 
                                         output_file=bins_seq_stats)
    if task_status:
        logger_ant.critical("Collecting sequences metrics FAILED!")
        logger_ant.critical("Terminating bins taxonomy identification and QC.")
        return 1, None
    logger_ant.info("Collecting sequences metrics COMPLETED.")
    logger_ant.debug(f"Collected metrics have been saved to '{bins_seq_stats}'.")
    
    logger_ant.debug(f"Masking bin_000.{bin_extension}...")
    try:
        os.rename(os.path.join(bins_dir, f"bin_000.{bin_extension}"),
                  os.path.join(bins_dir, "bin_000.MASKED"))
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)

    logger_ant.info("Identifying bins taxonomy...")
    gtdbtk_output_dir = os.path.join(work_dir, "gtdbtk")
    task_status = run_gtdbtk(bins_dir=bins_dir, 
                             output_dir=gtdbtk_output_dir, 
                             threads=threads,
                             bin_extension=bin_extension)
    if task_status:
        logger_ant.critical("Bins taxonomy indentification FAILED!")
        logger_ant.critical("Terminating bins taxonomy identification and QC.")
        return 1, None
    logger_ant.info("Bins taxonomy indentification COMPLETED.")

    logger_ant.info("Assessing bin completeness and contamination values...")
    if use_aa_seqs_for_qc:
        logger_ant.info("Using GTDB-Tk gene calls for bin quality metrics assessment...")
        bins_qc_input = os.path.join(work_dir, "gene_calls")
        os.makedirs(bins_qc_input)
        bins_qc_input_ext = "faa"
        logger_ant.info(f"Copying gene calls to the '{bins_qc_input}'...")
        move_gtdbtk_gene_calls(gtdbtk_output_dir=gtdbtk_output_dir, 
                               dest_dir=bins_qc_input,
                               ext=bins_qc_input_ext)
    else:
        logger_ant.info("Using bin nucleotide sequences for bin quality metrics assessment...")
        bins_qc_input = bins_dir
        bins_qc_input_ext = "fasta"
        
    task_status = bins_qc(bins_dir=bins_qc_input, 
                          output_dir=work_dir, 
                          threads=threads, 
                          bin_extension=bins_qc_input_ext)
    if task_status:
        logger_ant.critical("Assessment of bin completeness and contamination values FAILED!")
        logger_ant.critical("Terminating bins taxonomy identification and QC.")
        return 1, None
    logger_ant.info("Assessment of bin completeness and contamination values COMPLETED.")

    logger_ant.debug(f"Demasking bin_000.{bin_extension}...")
    try:
        os.rename(os.path.join(bins_dir, "bin_000.MASKED"),
                  os.path.join(bins_dir, f"bin_000.{bin_extension}"))
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)

    logger_ant.info("Combining data collected for bins...")
    task_status, MAGs = merge_tax_and_qc(gtdbtk_output_dir=gtdbtk_output_dir, 
                                         checkm_output_dir=os.path.join(work_dir, "checkm"), 
                                         checkm2_output_dir=os.path.join(work_dir, "checkm2"),  
                                         bins_stats_tsv=bins_seq_stats, 
                                         output_file=bins_metrics_tsv, 
                                         checkm_v=checkm_v)
    if task_status:
        logger_ant.warning("Combining data collected for bins FAILED!")
        return 1, None
    logger_ant.info("Combination of data collected for bins COMPLETED.")
    logger_ant.debug(f"Combined report have been saved to '{bins_metrics_tsv}'.")
    return 0, MAGs

##### Обработка результатов предсказания генов, генерация общих FASTA файлов и аннотации #####

def run_prodigal(genome_path, output_dir, 
                 mode="meta", output_format="gff", 
                 save_genes_nuc=True, save_genes_aa=True):
    '''
    Запускает функцию предсказание генов с использованием prodigal.
    
    Аргументы:
        genome_path (str)      путь к файлу с геномной последовательностью.
        output_dir (str)       путь к дирректории, куда будут сохраняться результаты работы prodigal.
        mode (str)             режим работы prodigal ["meta"].
        output_format (str)    формат генерируемой аннотации ["gff"].
        save_genes_nuc (bool)  сохранить нуклеотидные последовательности CDS? [True]
        save_genes_aa (bool)   сохранить аминокислотные последовательности CDS? [True] 
    '''
    env = "annotation"
    conda_run_cmd = f"conda run -n {env}"
    prodigal = "prodigal"
    fna = f"-d {os.path.join(output_dir, 'bin_000_genes.fna')}" if save_genes_nuc else ""
    faa = f"-a {os.path.join(output_dir, 'bin_000_genes.faa')}" if save_genes_aa else ""

    cmd = f"""{prodigal} -p {mode} 
              -f {output_format} {fna} {faa} 
              -i {genome_path} 
              -o {os.path.join(output_dir, f"bin_000_genes.{output_format}")}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def merge_fasta(fasta_paths, output_fasta, fix_prodigal_headers=True, 
                remove_description=False, add_to_description=[]):
    '''
    Объединяет несколько FASTA файлов в один.
    Если FASTA файлы были получены с использованием prodigal (последовательности CDS), 
    форматирует заголовки записей (опционально).
    Опционально, корректирует описание записей: удаляет полностью, либо добавляет доп. информацию.
    Собирает информацию об обнаруженных CDS и возвращает их число.
    
    Использует:
        - Bio.SeqIO
        
    Аргументы:
        fasta_paths (list)           список путей к FASTA файлам.
        output_fasta (str)           путь к файлу, куда будет записан результат.
        fix_prodigal_headers (bool)  переформатировать заголовки, сгенерированные prodigal? [True]
        remove_description (bool)    удалить описания записей? [False]
        add_to_description (list)    доп. информация для добавление в описание записей [[]].
                                     Желательно, список строк формата "key=value", но не обязательно.
    '''

    def fix_prodigal_fasta_headers(bio_seq_record):
        seq_id = ":".join(["_".join(bio_seq_record.id.split("_")[:-1]), bio_seq_record.id.split("_")[-1]])
        bio_seq_record.id = seq_id
        bio_seq_record.name = seq_id

        if remove_description:
            seq_description = ""
        else:
            seq_description = bio_seq_record.description.split(" # ")
            seq_description = [*seq_description[1:-1], *seq_description[-1].split(";")[1:]]
            seq_description[0] = f"start={seq_description[0]}"
            seq_description[1] = f"end={seq_description[1]}"
            seq_description[2] = f"strand={seq_description[2]}"
            seq_description = ";".join([f"ID={seq_id}", *seq_description, *add_to_description])
        bio_seq_record.description = seq_description
        return bio_seq_record

    n_cds = 0
    
    try:
        with open(output_fasta, 'a') as o_file:
            for fasta in fasta_paths:
                if os.path.exists(fasta):
                    with open(fasta) as i_file:
                        for seq in SeqIO.parse(i_file, "fasta"):
                            if fix_prodigal_headers:
                                seq = fix_prodigal_fasta_headers(seq)
                            SeqIO.write(seq, o_file, "fasta")
                            n_cds += 1
                else:
                    logger_ant.warning(f"File is not exists: '{fasta}'. Skip.")
                    continue
        return 0, n_cds
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1, None

def merge_gff(gff_paths, output_gff, fix_atr=True,
              add_to_atr=[], remove_comments=False):
    '''
    Объединяет несколько .gff файлов в один.
    Если .gff файлы были получены с использованием prodigal, форматирует строки атрибутов (опционально).
    Опционально, удаляет строки коментариев из финального файла; добавляет доп. информацию в атрибуты.
            
    Аргументы:
        gff_paths (list)        список путей к .gff файлам.
        output_gff (str)        путь к файлу, куда будет записан результат.
        fix_atr (bool)          переформатировать атрибуты, сгенерированные prodigal? [True]
        add_to_atr (list)       дополнительные атрибуты [[]].
                                Желательно, список строк формата "key=value", но не обязательно.
        remove_comments (bool)  удалить строки комметраиев? [False]
    '''
    try:
        with open(output_gff, 'a') as o_file:
            o_file.write("##gff-version\t3\n")
            for gff in gff_paths:
                if os.path.exists(gff):
                    with open(gff, 'r') as i_file:
                        for line in i_file:
                            if line[0] == "#":
                                if remove_comments or line[1] == "#":
                                    continue
                                o_file.write(line)
                            else:
                                if fix_atr:
                                    line = line.split()
                                    new_atr_line = f"ID={line[0]}:{line[-1].split("_", 1)[1]}"
                                    line = "\t".join([*line[:-1], new_atr_line])
                                if add_to_atr:
                                    line = f"{line}{";".join(add_to_atr)};"
                                o_file.write(f"{line}\n")
                else:
                    logger_ant.warning(f"File is not exists: '{fasta}'. Skip.")
                    continue
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def process_gene_calls(bins_dir, work_dir, output_dir,
                       bin_extension="fasta", save_cds_nuc_seqs=True):
    '''
    Реализует сведение результатов предсказания генов для отдельных бинов в общие файлы.
    1)  Выполняет предсказание генов для бина bin_000.
    2)  Проверяет существование директории с результатами предсказания генов для остальных бинов,
        выполненного ранее в процесе работы GTDB-Tk.
    3)  Объединяет и корректирует .gff, .faa и (опционально) .fna файлы отдельных бинов 
        в общие файлы и сохраняет их в отдельной директории.
    Собирает информацию об обнаруженных CDS и возвращает их число.

    Использует:
        - utils.py

    Аргументы:
        bins_dir (str)            путь к директории с последовательностями геномов.
        work_dir (str)            путь к рабочей дирректории, где хранятся результаты предсказания генов.
        output_dir (str)          путь к директории, в которой будет сохранены итоговые файлы.
        bin_extension (str)       расширение FASTA файлов ["fasta"].
        save_cds_nuc_seqs (bool)  сохранить нуклеотидные последовательности CDS? [True]
    '''
    logger_ant.debug("Bin 000 (unbinned contigs): predicting genes...")
    bin_000_gene_calls = os.path.join(work_dir, "bin_000_gene_calls")
    os.makedirs(bin_000_gene_calls)
    task_status = run_prodigal(genome_path=os.path.join(bins_dir, f"bin_000.{bin_extension}"), 
                               output_dir=bin_000_gene_calls, 
                               mode="meta", 
                               output_format="gff", 
                               save_genes_nuc=True, 
                               save_genes_aa=True)
    if task_status:
        logger_ant.critical("Bin 000 (unbinned contigs): genes prediction FAILED!")
        logger_ant.critical("Terminating gene calls processing.")
        return 1, None
    logger_ant.info("Bin 000 (unbinned contigs): genes prediction COMPLETED.")

    gc_dir = os.path.join(work_dir, "gtdbtk/identify/identify/intermediate_results/marker_genes")
    if not os.path.exists(gc_dir):
        logger_ant.critical("Unable to find directory with gene calls!")
        logger_ant.critical("Terminating gene calls processing.")
        return 1, None
    logger_ant.info(f"Directory '{gc_dir}' found.")
    gc_file_list = [os.path.join(gc_dir, mg_bin, f"{mg_bin}_protein") for mg_bin in os.listdir(gc_dir)]
    gc_file_list.append(os.path.join(bin_000_gene_calls, "bin_000_genes"))

    got_warnings = False
    if len(gc_file_list) != len(os.listdir(bins_dir)):
        logger_ant.warning("The number of gene calling outputs found does not match the number of bins!")
        logger_ant.warning(f"Found - {len(gc_file_list)}; expected - {len(os.listdir(bins_dir))}.")
        got_warnings = True

    logger_ant.info("Merging .gff files...")
    mg_assembly_annot_file = os.path.join(output_dir, "cds_annot.gff")
    task_status = merge_gff(gff_paths=[f"{path}.gff" for path in gc_file_list], 
                            output_gff=mg_assembly_annot_file, 
                            fix_atr=True,
                            add_to_atr=[], 
                            remove_comments=False)
    if task_status:
        logger_ant.critical("Merging .gff files FAILED!")
        logger_ant.critical("Terminating gene calls processing.")
        return 1, None
    logger_ant.info("Merging .gff files COMPLETED.")
    logger_ant.info(f"Path to the merged .gff file: {mg_assembly_annot_file}.")

    logger_ant.info("Merging .fasta files with amino acid sequences of CDS...")
    cds_faa_file = os.path.join(output_dir, "cds_aa_seqs.faa")
    task_status, n_cds = merge_fasta(fasta_paths=[f"{path}.faa" for path in gc_file_list], 
                                     output_fasta=cds_faa_file, 
                                     fix_prodigal_headers=True, 
                                     remove_description=False, 
                                     add_to_description=[])
    if task_status:
        logger_ant.critical("Merging .fasta files with amino acid sequences of CDS FAILED!")
        logger_ant.critical("Terminating gene calls processing.")
        return 1, None
    logger_ant.info("Merging .fasta files with amino acid sequences of CDS COMPLETED.")
    logger_ant.info(f"Path to the merged .faa file: {cds_faa_file}.")

    if save_cds_nuc_seqs:
        logger_ant.info("Merging .fasta files with nucleotide sequences of CDS...")
        cds_fna_file = os.path.join(output_dir, "cds_nuc_seqs.fna")
        task_status, _ = merge_fasta(fasta_paths=[f"{path}.fna" for path in gc_file_list], 
                                     output_fasta=cds_fna_file, 
                                     fix_prodigal_headers=True, 
                                     remove_description=False, 
                                     add_to_description=[])
        if task_status:
            logger_ant.warning("Merging .fasta files with nucleotide sequences of CDS FAILED!")
            got_warnings = True
        else:
            logger_ant.info("Merging .fasta files with nucleotide sequences of CDS COMPLETED.")
            logger_ant.info(f"Path to the merged .fna file: {cds_fna_file}.")

    if got_warnings:
        return 2, n_cds
    return 0, n_cds

##### Аннотация CDS c MMseqs2 #####

def run_mmseqs2_createdb(fasta_file, seqs_db, threads):
    '''
    Запускает функцию createdb пакета MMseqs2. Создает БД последовательностей.

    Аргументы:
        fasta_file (str)  путь к файлу с последовательностями.
        seqs_db (str)     путь к дирректории, куда будет записана БД, в формате /path/to/dir/db_basename.
        threads (int)     число потоков.
    '''
    env = "mmseqs2"
    conda_run_cmd = f"conda run -n {env}"
    mmseqs = "mmseqs"

    cmd = f"{mmseqs} createdb {fasta_file} {seqs_db} --threads {threads}"
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_mmseqs2_createsubdb(subset_info, input_db, output_db):
    '''
    Запускает функцию createsubdb пакета MMseqs2. Создает субБД последовательностей.

    Аргументы:
        subset_info (str)  путь к БД, на основе которой будет сформирована подвыборка последобательностей.
                           Непример, результаты кластеризации последовательностей.
                           Формат: /path/to/dir/db_basename.
        input_db (str)     путь к исходной БД. Формат: /path/to/dir/db_basename.
        output_db (str)    путь к исходной субБД. Формат: /path/to/dir/db_basename.
    '''
    env = "mmseqs2"
    conda_run_cmd = f"conda run -n {env}"
    mmseqs = "mmseqs"

    cmd = f"{mmseqs} createsubdb {subset_info} {input_db} {output_db}"
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_mmseqs2_cluster(input_db, cluster_db, tmp_dir, threads, 
                        seq_ident="90%"):
    '''
    Запускает функцию createsubdb пакета MMseqs2. Создает субБД последовательностей.

    Аргументы:
        input_db (str)    путь к исходной БД. Формат: /path/to/dir/db_basename.
        cluster_db (str)  путь к БД, куда будут записаны результаты кластеризации. 
                          Формат: /path/to/dir/db_basename.
        tmp_dir (str)     путь к директории для хранения временных файлов.
        threads (int)     число потоков.
        seq_ident (str)   порог идентичности последовательностей ["90%"].
                          Принитмает значения "90%" или "50%".
    '''
    env = "mmseqs2"
    conda_run_cmd = f"conda run -n {env}"
    mmseqs = "mmseqs"

    ident_settings = {"90%" : [0.9, 0.95],
                      "50%" : [0.5, 0.8]}

    cmd = f"""{mmseqs} cluster {input_db} {cluster_db} {tmp_dir} 
              --min-seq-id {ident_settings[seq_ident][0]} 
              -c {ident_settings[seq_ident][1]} 
              --cov-mode 1 
              --threads {threads}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_mmseqs2_createtsv(query_db, target_db, result_db, output_tsv, threads):
    '''
    Запускает функцию createtsv пакета MMseqs2. 
    Генерирует .tsv файл, отражающий взаимосвязи между записямиразными БД.
    Например, генерирует .tsv файл, отражающий взаимосвязь между исходными последовательностями (query_db) и
    последовательностями в референсной БД (target_db) на основе результатов аннотации (result_db).

    Аргументы:
        query_db (str)    путь к БД. Формат: /path/to/dir/db_basename.
        target_db (str)   путь к БД. Формат: /path/to/dir/db_basename.
        result_db (str)   путь к БД. Формат: /path/to/dir/db_basename.
        output_tsv (str)  путь к генерируемому .tsv файлу.
        threads (int)     число потоков.
    '''
    env = "mmseqs2"
    conda_run_cmd = f"conda run -n {env}"
    mmseqs = "mmseqs"

    cmd = f"{mmseqs} createtsv {query_db} {target_db} {result_db} {output_tsv} --threads {threads}"
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        logger_ant.info("You can try to create this table manualy using 'createtsv' MMseqs2 function.")
        return 1

def run_mmseqs2_convert2fasta(input_db, output_fasta):
    '''
    Запускает функцию convert2fasta пакета MMseqs2. 
    Генерирует FASTA файл с последовательностями на сонове БД.

    Аргументы:
        input_db (str)      путь к БД. Формат: /path/to/dir/db_basename.
        output_fasta (str)  путь к генерируемому FASTA файлу.
    '''
    env = "mmseqs2"
    conda_run_cmd = f"conda run -n {env}"
    mmseqs = "mmseqs"

    cmd = f"{mmseqs} convert2fasta {input_db} {output_fasta}"
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred.", exc_info=True)
        logger_ant.info("You can try to create this file manualy using 'convert2fasta' MMseqs2 function.")
        return 1

def run_mmseqs2_search(query_db, ref_db, result_db, tmp_dir, threads):
    '''
    Запускает функцию search пакета MMseqs2. 
    Осужествляет поиск последовательностей из одной БД (query_db) в другой (ref_db).

    Аргументы:
        query_db (str)   путь к исходной БД. Формат: /path/to/dir/db_basename.
        ref_db (str)     путь к референсной БД. Формат: /path/to/dir/db_basename.
        result_db (str)  путь к БД для записи результатов поиска. Формат: /path/to/dir/db_basename.
        tmp_dir (str)    путь к директории для хранения временных файлов.
        threads (int)    число потоков.
    '''
    env = "mmseqs2"
    conda_run_cmd = f"conda run -n {env}"
    mmseqs = "mmseqs"

    cmd = f"{mmseqs} search {query_db} {ref_db} {result_db} {tmp_dir} --threads {threads}"
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def generate_annotation_table(input_db, ref_db, annotation_db, annotation_tsv, tmp_dir, threads, 
                              only_first_occurrences=True):
    '''
    Формирует таблицу аннотации на на основе .tsv файла, генерируемого run_mmseqs2_createtsv для
    результатов аннотации CDS против референсной ДБ.
    Опционально, фильтрует ее, отсавляя только первые (наиболее достоверные) вхождения для каждой CDS.

    Аргументы:
        input_db (str)                путь к исходной БД. Формат: /path/to/dir/db_basename.
        ref_db (str)                  путь к референсной БД. Формат: /path/to/dir/db_basename.
        annotation_db (str)           путь к БД с результатами аннотации. Формат: /path/to/dir/db_basename.
        annotation_tsv (str)          путь к .tsv файлу, в который будет записана таблица аннотации.
        tmp_dir (str)                 путь к директории для хранения временных файлов.
        threads (int)                 число потоков.
        only_first_occurrences (bool) оставить только первые вхождения для каждой CDS? [True]
    '''
    annotation_tsv_full = os.path.join(tmp_dir, f"{annotation_tsv.split("/")[-1]}.full") if only_first_occurrences else annotation_tsv
    task_status = run_mmseqs2_createtsv(query_db=input_db, 
                                        target_db=ref_db, 
                                        result_db=annotation_db, 
                                        output_tsv=annotation_tsv_full, 
                                        threads=threads)
    if task_status:
        logger_ant.warning(f"Unable to create .tsv file with annotation from the {annotation_db.split("/")[-1]}.")
        return 1

    if only_first_occurrences:
        logger_ant.info("Removing not the first occurrences...")
        cmd = f"awk -F'\t' '!seen[$1]++' {annotation_tsv_full} > {annotation_tsv}"
        logger_ant.debug(f"CMD: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            logger_ant.info("Not the first occurrences removed successfully.")
            return 0
        except subprocess.CalledProcessError as error:
            logger_ant.error("An error occurred!", exc_info=True)
            return 1

def annotate_genes(cds_fasta_file, ref_db, work_dir, output_dir, threads, 
                   search_rep_seqs=True, seq_ident="90%", only_first_occurrences=True):
    '''
    Реализует пайплайн аннотации белковых последовательностей MMseqs2.
    1)  (Опционально) кластеризует аминокислотные последовательности для выявления репрезентативных.
        Сохраняет .tsv файл с отношениями seq - rep_seq в output_dir.
    2)  Выполняет поиск последовательностей в референсной БД (UniRef90).
    3)  Генерирует таблицу аннотации в .tsv формате и сохраняет ее в output_dir.

    Аргументы:
        cds_fasta_file (str)           путь к FASTA файлу с последовательностями CDS.
        ref_db (str)                   путь к референсной БД (UniRef90).
        work_dir (str)                 путь к директории, в которой будут сохранены результаты работы MMseqs2.
        output_dir (str)               путь к директории, в которой будут сохранены итоговые файлы.
        threads (int)                  число потоков.
        search_rep_seqs (bool)         работать в репрезентативные последовательностями? [True]
                                       Если нет, то этап кластеризации последовательностей будет пропущен и 
                                       в работу будут взяты все предоставленные последовательности.
        seq_ident (str)                порог идентичности последовательностей ["90%"].
                                       Принитмает значения "90%" или "50%".
        only_first_occurrences (bool)  оставить только первые вхождения для каждой CDS? [True]
                                       Если нет, то таблица аннотации будет включать все обнаруженные 
                                       достоверные совпадения (до 300 на последовательность).
    '''
    mmseqs2_dir = os.path.join(work_dir, "mmseqs2")
    os.makedirs(mmseqs2_dir)
 
    logger_ant.info("Creating CSD sequeces database...")
    seqs_db = os.path.join(mmseqs2_dir, "seqs_db")
    task_status = run_mmseqs2_createdb(fasta_file=cds_fasta_file, 
                                       seqs_db=seqs_db, 
                                       threads=threads)
    if task_status:
        logger_ant.critical("CSD sequeces database creation FAILED!")
        logger_ant.critical("Terminating CDS annotation.")
        return 1

    if search_rep_seqs:
        logger_ant.info("Working with representative sequences.")
        logger_ant.info("Clustering CSD sequences...")
        cluster_db = os.path.join(mmseqs2_dir, "cluster_db")
        task_status = run_mmseqs2_cluster(input_db=seqs_db, 
                                          cluster_db=cluster_db, 
                                          tmp_dir=os.path.join(mmseqs2_dir, "tmp"),
                                          threads=threads,
                                          seq_ident=seq_ident)
        if task_status:
            logger_ant.critical("CSD sequences clusterization FAILED!")
            logger_ant.critical("Terminating CDS annotation.")
            return 1

        got_warnings = False
        geneid2rep_tsv = os.path.join(output_dir, "geneid2rep.tsv")
        task_status = run_mmseqs2_createtsv(query_db=seqs_db, 
                                            target_db=seqs_db, 
                                            result_db=cluster_db, 
                                            output_tsv=geneid2rep_tsv,
                                            threads=threads)
        if task_status:
            logger_ant.warning("Unable to save gene ID to representative sequence relations data.")
            got_warnings = True
        else:
            logger_ant.info(f"Gene ID to representative sequence relations data was written to the '{geneid2rep_tsv}'.")

        logger_ant.info("Creating a database with representative cluster sequences...")
        rep_seqs_db = os.path.join(mmseqs2_dir, "rep_seqs_db")
        task_status = run_mmseqs2_createsubdb(subset_info=cluster_db, 
                                              input_db=seqs_db, 
                                              output_db=rep_seqs_db)
        if task_status:
            logger_ant.critical("Database with representative cluster sequences creation FAILED!")
            logger_ant.critical("Terminating CDS annotation.")
            return 1
        seqs_db = rep_seqs_db
        logger_ant.info("Database with representative cluster sequences creation COMPLETED.")
    else:
        logger_ant.info("Working with all sequences.")
    
    logger_ant.info("Annotation of CDS using UniRef90 database...")
    logger_ant.info("Searching sequences...")
    annotation_db = os.path.join(mmseqs2_dir, "UniRef90_annot_db")
    task_status = run_mmseqs2_search(query_db=seqs_db, 
                                     ref_db=ref_db, 
                                     result_db=annotation_db, 
                                     tmp_dir=os.path.join(mmseqs2_dir, "tmp"), 
                                     threads=threads)
    if task_status:
        logger_ant.warning("Sequences search FAILED!")
        return 2
    logger_ant.info("Sequences search COMPLETED.")
        
    logger_ant.info("Creating annotation table...")
    annotation_tsv = os.path.join(output_dir, "cds2UniRef90_annot.tsv")
    task_status = generate_annotation_table(input_db=seqs_db, 
                                            ref_db=ref_db, 
                                            annotation_db=annotation_db, 
                                            annotation_tsv=annotation_tsv, 
                                            tmp_dir=os.path.join(mmseqs2_dir, "tmp"),
                                            threads=threads, 
                                            only_first_occurrences=only_first_occurrences)
    if task_status:
        logger_ant.warning("The creation of the annotation table FAILED!")
        return 2
    logger_ant.info("The creation of the annotation table COMPLETED.")
    logger_ant.info(f"Annotation table path: {annotation_tsv}.")

    if got_warnings:
        return 2
    return 0

##### Генерация ридкаунтов для CDS #####

def run_bwa_index(assembly_fasta, output_dir):
    '''
    Запускает функцию index пакета bwa. Создает индекс последовательности.

    Аргументы:
        assembly_fasta (str)  путь к FASTA файлу с нуклеотидной последовательсностью.
        output_dir (str)      путь к дирректории, куда будет сохранен индекс.
    '''
    env = "annotation"
    conda_run_cmd = f"conda run -n {env}"
    bwa = "bwa"

    prefix = os.path.join(output_dir, "index")
    cmd = f"{bwa} index -p {prefix} {assembly_fasta}"
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_bwa_mem(R1, R2, index_path, output_file, threads):
    '''
    Запускает функцию mem пакета bwa. Производит выравнивание ридов.
    Если R1 и R2 асинхронны, пытается исправить это путем вызова функции run_bbmap_repair
    модуля utils, после чего делает повторную попытку выравнивания.

    Использует:
        - utils.py

    Аргументы:
        R1 (str)           путь к .fastq файлу предобработанных R1 ридов.
        R2 (str)           путь к .fastq файлу предобработанных R2 ридов.
        index_path (str)   путь к индексу в формате /path/to/dir/index.
        output_file (str)  путь к файлу выравнивания в формате .sam.
        threads (int)      число потоков.
    '''
    
    env = "annotation"
    conda_run_cmd = f"conda run -n {env}"
    bwa = "bwa"

    cmd = f"bwa mem -t {threads} -o {output_file} {index_path} {R1} {R2}"
    logger_ant.debug(f"CMD: {" ".join(cmd.split())}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", 
                       shell=True, capture_output=True, text=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        process_output = error.stderr
        logger_ant.error("An error occurred!", exc_info=True)

    error_message = "[mem_sam_pe] paired reads have different names:"
    if error_message in process_output:
        logger_ant.warning("Looks like the .fastq files provided are asynchronous.")
        logger_ant.info("Trying to fix it...")
    else:
        return 1
        
    R1_repaired = R1.replace(".fastq.gz", "_repaired.fastq.gz")
    R2_repaired = R2.replace(".fastq.gz", "_repaired.fastq.gz")
    task_status = utils.run_bbmap_repair(r1_input=R1, 
                                         r2_input=R2, 
                                         r1_output=R1_repaired, 
                                         r2_output=R2_repaired)
    if task_status:
        logger_ant.error("Unable to fix the provided .fastq files!")
        return 1

def run_coverm(bam_file, mode, output_file, threads,
               metrics=["count"], bins_dir=None, ext="fasta",
               exclude_supp=True):
    '''
    Запускает функции contig и genome пакета CoverM. Проводит подсчет ридов и метрик покрытия
    контигов и бинов на основе выравниваний.

    Аргументы:
        assembly_fasta (str)  путь к FASTA файлу с нуклеотидной последовательсностью.
        output_dir (str)      путь к дирректории, куда будет сохранен индекс.
    '''
    env = "annotation"
    conda_run_cmd = f"conda run -n {env}"
    coverm = "coverm"

    if mode == "genome" and not bins_dir:
        logger_ant.error("In 'genome' mode the path to directory with bins must be specified!")
        return 1

    metrics_avail = {"contig" : ["mean", "trimmed_mean", "coverage_histogram",
                                 "covered_fraction", "covered_bases", "variance", 
                                 "length", "count", "metabat", "reads_per_base", 
                                 "rpkm", "tpm"],
                     "genome" : ["relative_abundance", "mean", "trimmed_mean",
                                 "coverage_histogram", "covered_fraction", "covered_bases",
                                 "variance", "length", "count", "reads_per_base", 
                                 "rpkm", "tpm"]}
    metrics_accept = [m for m in metrics if m in metrics_avail[mode]]
    metrics_disc = [m for m in metrics if m not in metrics_avail[mode]]
    if metrics_disc:
        logger_ant.warning(f"Metrics {" ".join(metrics_disc)} are not accepted by the CoverM and will be ignored.")

    if mode == "contig":
        cmd = f"""{coverm} {mode} -t {threads} 
                  -m {" ".join(metrics_accept)} 
                  -b {bam_file} 
                  -o {output_file}"""
    elif mode == "genome":
        ex_sup_reads = "--exclude-supplementary" if exclude_supp else ""
        cmd = f"""{coverm} {mode} -t {threads} 
                  -m {" ".join(metrics_accept)} 
                  {ex_sup_reads} 
                  --min-covered-fraction 0
                  -b {bam_file}
                  -d {bins_dir}
                  -x {ext}
                  -o {output_file}"""
    else:
        logger_ant.error("Unexpected CoverM mode!")
        return 1

    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def run_featurecounts(alignment, gff_path, output_file, threads,
                      feature_type="CDS", feature_id="ID", min_aln_q=20):
    '''
    Запускает функцию featureCounts пакета Subread. Производит подсчет ридов.
    Принимает во внимание только первичные выравнивания. 
    Парные риды считает по отдельности.
    Для мультимэперов и мультиоверлэперов присваивает дробный каунт 
    в соответствии с числом охватываемых фич. 

    Аргументы:
        alignment (str)     путь к файлу выравнивания.
        gff_path (str)      путь к аннотации в .gff формате.
        output_file (str)   путь к файлу, в который будут записаны ридкаунты.
        threads (int)       число потоков.
        feature_type (str)  тип фичи, по которой будет производится подсчет ридов ["CDS"].
        feature_id (str)    идентификатор фичи ["ID"].
        min_aln_q (int)     минимальное допустимое качество выравнивания [20].
    '''
    
    env = "annotation"
    conda_run_cmd = f"conda run -n {env}"
    feature_counts = "featureCounts"

    cmd = f"""{feature_counts} -T {threads} 
              -p --primary -O -M --fraction -Q {min_aln_q} 
              -t {feature_type} -g {feature_id} 
              -a {gff_path} 
              -o {output_file} 
              {alignment}"""
    cmd = " ".join(cmd.split())
    logger_ant.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_ant.error("An error occurred!", exc_info=True)
        return 1

def get_readcounts_for_features(assembly_fasta, gff_file, reads_dir, bins_dir, work_dir, output_dir, threads,
                                bin_extension="fasta", feature_type="CDS", feature_id="ID", min_aln_q=20):
    '''
    Реализует все этапы процесса подсчета ридов, выравненных на CDS, обнаруженных.
    1)  Выполняет и процессирует выравнивания ридов. Отчет со статистикой выравнивания
        в tsv. формате сохраняет в output_dir.
    2)  Выполняет подсчет ридов пересекающихся с отдельными контигами, бинами и CDS 
        и сохраняет файл с ридкаунтами и отчет о подсчете в output_dir.
    В процессе собирается статистика выравнивания и подсчета ридов, 
    которая возвращается в виде словаря со следующей структурой:
    aln_stats = {
        "readsTotal" : int,
        "mappedReads" : int,
        "mappedReadsPrcnt" : float,
        "assignedReads" : int,
        "assignedReadsPrcnt" : float
    }

    Использует:
        - utils.py

    Аргументы:
        assembly_fasta (str)  путь к FASTA файлу с последовательностью метагеномной сборки.
        gff_file (str)        путь к .gff файлу аннотации метагеномной сборки.
        reads_dir (str)       путь к директории, в которой хранятся предобработанные риды.
        bins_dir (str)        путь к директории с FASTA файлами с последовательностями бинов.
        work_dir (str)        путь к директории, в которой будут сохранены промежуточные результаты.
        output_dir (str)      путь к директории, в которой будут сохранены итоговые файлы.
        threads (int)         число потоков.
        bin_extension (str)   расширение FASTA файлов с последовательностями бинов ["fasta"].
        feature_type (str)    тип фичи, по которой будет производится подсчет ридов ["CDS"].
        feature_id (str)      идентификатор фичи ["ID"].
        min_aln_q (int)       минимальное допустимое качество выравнивания [20].
    '''
    required_files = [assembly_fasta, gff_file,
                      os.path.join(reads_dir, "R1_paired.fastq.gz"),
                      os.path.join(reads_dir, "R2_paired.fastq.gz")]
    for required_file in required_files:
        if os.path.exists(required_file):
           logger_ant.info(f"File '{required_file}' found.")
        else:
            logger_ant.critical(f"The required '{required_file}' does not exist!")
            logger_ant.critical("Terminating reads counting.")
            return 1, None
            
    rc_dir = os.path.join(work_dir, "read_counting")
    os.makedirs(rc_dir, exist_ok=True)

    logger_ant.info("Building assembly index...")
    task_status = utils.run_bowtie2_build(fasta_file=assembly_fasta, 
                                          output_dir=rc_dir, 
                                          threads=threads)
    if task_status:
        logger_ant.critical("Building assembly index FAILED!")
        logger_ant.critical("Terminating reads counting.")
        return 1, None
    logger_ant.info("Assembly index building COMPLETED.")

    logger_ant.info("Aligning...")
    task_status = utils.run_bowtie2(reads_dir=reads_dir, 
                                    index_path=os.path.join(rc_dir, "index"), 
                                    output_file=os.path.join(rc_dir, "alignment.sam"), 
                                    threads=threads,
                                    use_unpaired=False)
    if task_status:
        logger_ant.critical("Alignment FAILED!")
        logger_ant.critical("Terminating reads counting.")
        return 1, None

    logger_ant.info("Processing alignments...")
    logger_ant.info("Converting .sam to .bam...")
    task_status = utils.convert_sam_to_bam(sam_file=os.path.join(rc_dir, "alignment.sam"), 
                                           bam_file=os.path.join(rc_dir, "alignment_unsorted.bam"), 
                                           threads=threads, 
                                           exclude_unmapped=True)
    if task_status:
        logger_ant.critical("Conversion .sam to .bam FAILED!")
        logger_ant.critical("Terminating reads counting.")
        return 1, None
    logger_ant.info("Sorting .bam...")
    task_status = utils.sort_bam(unsorted_bam_file=os.path.join(rc_dir, "alignment_unsorted.bam"), 
                                 sorted_bam_file=os.path.join(rc_dir, "alignment.bam"), 
                                 threads=threads, 
                                 by_name=False)
    if task_status:
        logger_ant.critical("Sorting FAILED!")
        logger_ant.critical("Terminating reads counting.")
        return 1, None
    logger_ant.info("Aligning COMPLETED.")

    logger_ant.info("Getting flag statistics for the alignment...")
    flagstat_report = os.path.join(output_dir, "reads2assembly_alignment_stats.tsv")
    task_status = utils.get_alignment_flagstat(aln_file=os.path.join(rc_dir, "alignment.sam"), 
                                               report_file=flagstat_report,
                                               report_format="tsv",
                                               threads=threads)
    if task_status:
        logger_ant.warning("Can't get flag statistics!")
    else:
        logger_ant.info(f"Flag statistics was saved to '{flagstat_report}'.")
    logger_ant.info("Alignments processing COMPLETED.")

    logger_ant.info("Per contig reads counting...")
    metrics = ["length", "count", "mean", "trimmed_mean", "covered_fraction", 
               "covered_bases", "reads_per_base", "rpkm", "tpm"]
    task_status = run_coverm(bam_file=os.path.join(rc_dir, "alignment.bam"), 
                             mode="contig", 
                             output_file=os.path.join(output_dir, "reads_per_contig.tsv"), 
                             threads=threads,
                             metrics=metrics)
    if task_status:
        logger_ant.critical("Per contig reads counting FAILED.")
        logger_ant.critical("Terminating reads counting.")
        return 1, None

    logger_ant.info("Per bin reads counting...")
    metrics = ["length", "count", "relative_abundance", "mean", "trimmed_mean",
               "covered_fraction", "covered_bases", "reads_per_base", "rpkm", "tpm"]
    task_status = run_coverm(bam_file=os.path.join(rc_dir, "alignment.bam"), 
                             mode="genome", 
                             output_file=os.path.join(output_dir, "reads_per_bin.tsv"), 
                             threads=threads,
                             metrics=metrics,
                             bins_dir=bins_dir, 
                             ext=bin_extension,
                             exclude_supp=True)
    if task_status:
        logger_ant.critical("Per bin reads counting FAILED.")
        logger_ant.critical("Terminating reads counting.")
        return 1, None

    logger_ant.info("Per CDS reads counting...")
    task_status = run_featurecounts(alignment=os.path.join(rc_dir, "alignment.bam"), 
                                    gff_path=gff_file, 
                                    output_file=os.path.join(output_dir, "reads_per_cds.tsv"),
                                    threads=threads,
                                    feature_type=feature_type, 
                                    feature_id=feature_id, 
                                    min_aln_q=min_aln_q)
    if task_status:
        logger_ant.critical("Counting reads FAILED.")
        logger_ant.critical("Terminating reads counting.")
        return 1, None
    featurecounts_report = os.path.join(output_dir, "featurecounts_report.tsv")
    os.rename(os.path.join(output_dir, "reads_per_cds.tsv.summary"), featurecounts_report)
    logger_ant.info(f"Read counting statistics was saved to '{featurecounts_report}'.")

    logger_ant.info("Collecting alignment metrics...")
    aln_stats = {}
    try:
        flagstats = pd.read_csv(flagstat_report, sep="\t", header=None)
        reads_total = int(flagstats.iloc[0, 0])
        reads_map = int(flagstats.loc[flagstats[2] == "mapped", 0].iloc[0])
        aln_stats["readsTotal"] = reads_total
        aln_stats["mappedReads"] = reads_map
        aln_stats["mappedReadsPrcnt"] = round(reads_map / reads_total * 100, 2)
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)
        logger.warning("Unable to collect alignmnet metrics!")
    try:
        with open(featurecounts_report, "r") as report:
            for line in report:
                if "Assigned" in line:
                    reads_assign = int(line.split()[1])
                    break
        aln_stats["assignedReads"] = reads_assign
        aln_stats["assignedReadsPrcnt"] = round(reads_assign / reads_total * 100, 2)
    except Exception as error:
        logger_ant.error("An error occurred!", exc_info=True)
        logger.warning("Unable to collect read counting metrics!")
    
    logger_ant.info("Reads counting COMPLETED.")
    return 0, aln_stats

##### Главный скрипт модуля аннотации #####

def annotate_mg(sample, threads, ref_db,
                use_aa_seqs_for_bins_qc=True, 
                checkm_v_for_mags_q_asses=1,
                save_cds_nuc_seqs=True,
                use_repr_seqs_for_cds_annot=True,
                clean_up=True):
    '''
    Реализует полный процесс аннотации метагеномной сборки: 
    0)  Проверяет существование требуемых файлов.
    1)  Проводит таксономическую идентификацию и контроль качества бинов.
        Результат сохранятеся в виде сводной таблицы в файл bin_metrics.tsv в директории assembly.
    2)  Осуществляет постпроцессинг файлов с аннотациями и последовательностями генов,
        полученными на предыдущем этапе. Полученные файлы с аннотацией всей метагеномной сборки (cds_annot.gff),
        аминокислотными (cds_aa_seqs.faa) и нуклеотидными последовательностями (cds_nuc_seqs.fna, опционально)
        CDS сохраняются в директории annotations.
    3)  Аннотирует обнаруженные CDS против БД UniRef90 с иcпользованием пайплайна MMseqs2.
        Полученная таблица аннотации (cds2UniRef90_annot.tsv) сохраняется в директории annotations.
        Если аннотация выполнялась только для рперезентативных последовательностей, в директории annotations
        также сохраняется файл с отношениями seq - rep_seq (geneid2rep.tsv).
    4)  Проводит подсчет ридов, относящихся к отдельным контигам, бинам и CDS. 
        Файлы с ридкаунтами, а также репорты со статистикой выравнивания и подсчета ридов 
        сохраняются в директории readcounts.
    В процессе работы собирается метрики: количество MAGs разных категорий качества, число обнаруженных CDS,
    статистики выравнивания и подсчета ридов.

    В результате работы в директории сборки образца остаются:
    sample_dir
    |----- assembly.log
    |----- assembly
    |      |----- contigs.fasta
    |      |----- assembly_graph.fastg
    |      |----- contigs.paths
    |      |----- scaffolds.paths
    |      |----- bin_metrics.tsv
    |      |----- contig2bin.tsv
    |      |----- bins
    |             |----- bin_000.fasta
    |             |----- bin_0.fasta
    |             |----- ...
    |----- reports
    |      |----- spades.log
    |      |----- fastp_report.html
    |      |----- fastp_report.json
    |----- annotations
    |      |----- cds_annot.gff
    |      |----- cds_aa_seqs.faa
    |      |----- cds_nuc_seqs.fna (опционально)
    |      |----- geneid2rep.tsv (опционально)
    |      |----- cds2UniRef90_annot.tsv
    |----- readcounts
           |----- reads_per_contig.tsv
           |----- reads_per_bin.tsv
           |----- reads_per_cds.tsv
           |----- reads2assembly_alignment_stats.tsv
           |----- featurecounts_report.tsv
    
    Использует:
        - sample_cl.py
        - utils.py
        
    Аргументы:
        sample                              экземпляр класса AnnotatedSample, созданный для конкретного образца,
                                            определяющий структуру содержания директории сборки образца и 
                                            хранящий информацию о статусе аннотации, файлах с важными результатами,
                                            а также собираемые метрики.
        ref_db (str)                        Путь к рефернсной БД (принимаются готовые БД под MMseqs2),
                                            на основе которой будут аннотироваться CDS.
        threads (int)                       число потоков.
        use_aa_seqs_for_bins_qc (bool)      Использовать аминокислотные последовательности CDS для QC бинов? 
                                            Если нет, то QC бинов будет осуществляться на полных 
                                            нуклеотидных последовательностях бинов.
        checkm_v_for_mags_q_asses (bool)    версия CheckM, чьи метрики завершенности и контаминации использовать
                                            для присвоения категории качества MAGs [1]. CheckM - 1; CheckM2 - 2.
        save_cds_nuc_seqs (bool)            сохранить нуклеотидные последовательности CDS? [True]
        use_repr_seqs_for_cds_annot (bool)  работать в репрезентативные последовательностями? [True]
                                            Если нет, то этап кластеризации последовательностей будет пропущен и 
                                            в работу будут взяты все предоставленные последовательности.
        clean_up (bool)                     удалить все промежуточные файлы? [True] 

    Возвращает:
        - статус сборки образца ["Annotated" | "Annotation failed"]
        - словарь с информацией о статусах подпроцессов, файлах с результатами и метриками.
    '''
    logger_ant.debug(f"Working in '{os.getcwd()}'.")
    
    logger_ant.info("Verifying inputs...")
    try:
        sample.check_inputs()
    except FileNotFoundError as error:
        logger_ant.error("Inputs verification FAILED!", exc_info=True)
        logger_ant.critical("Terminating sample processing.")
        sample.set_process_status("inputsCheck", "FAILED")
        sample.status = "Annotation failed"
        return sample
    sample.set_process_status("inputsCheck", "COMPLETED")
    logger_ant.info("Inputs verification COMPLETED.")
    
    os.makedirs(sample.swd)
    logger_ant.info("Starting bins taxonomy identification and QC...")
    task_status, MAGs = bins_tax_ident_and_qc(bins_dir=sample.bins, 
                                              work_dir=sample.swd, 
                                              bins_metrics_tsv=sample.outputs["binsMetrics"],
                                              threads=threads, 
                                              bin_extension="fasta", 
                                              use_aa_seqs_for_qc=use_aa_seqs_for_bins_qc, 
                                              checkm_v=checkm_v_for_mags_q_asses)
    if task_status:
        logger_ant.critical("Bins taxonomy identification and QC FAILED!")
        logger_ant.critical("Terminating sample processing.")
        sample.set_process_status("taxonomyIdent&QC", "FAILED")
        sample.status = "Annotation failed"
        return sample
    logger_ant.info("Bins taxonomy identification and QC COMPLETED.")
    sample.memoize_metric("MAGs", MAGs)
    sample.set_process_status("taxonomyIdent&QC", "COMPLETED")

    os.makedirs(sample.annots)
    logger_ant.info("Starting processing gene calls...")
    task_status, n_cds = process_gene_calls(bins_dir=sample.bins, 
                                            work_dir=sample.swd, 
                                            output_dir=sample.annots,
                                            bin_extension="fasta", 
                                            save_cds_nuc_seqs=save_cds_nuc_seqs)
    if task_status:
        status = "FAILED" if task_status == 1 else "FINISHED WITH WARNINGS"
        logger_ant.critical(f"Gene calls processing {status}!")
        logger_ant.critical("Terminating sample processing.")
        sample.set_process_status("geneCallsProcessing", status)
        sample.status = "Annotation failed"
        return sample
    logger_ant.info("Gene calls processing COMPLETED.")
    sample.memoize_metric("CDS", n_cds)
    sample.set_process_status("geneCallsProcessing", "COMPLETED")

    logger_ant.info("Starting CDS annotation...")
    task_status = annotate_genes(cds_fasta_file=sample.outputs["CDSaaFasta"], 
                                 ref_db=ref_db, 
                                 work_dir=sample.swd, 
                                 output_dir=sample.annots, 
                                 threads=threads, 
                                 search_rep_seqs=use_repr_seqs_for_cds_annot, 
                                 seq_ident="90%", 
                                 only_first_occurrences=True)
    if task_status:
        status = "FAILED" if task_status == 1 else "FINISHED WITH WARNINGS"
        logger_ant.critical(f"Gene calls processing {status}!")
        logger_ant.critical("Terminating sample processing.")
        sample.set_process_status("annotating", status)
        sample.status = "Annotation failed"
        return sample
    logger_ant.info("CDS annotation COMPLETED.")
    sample.set_process_status("annotating", "COMPLETED")

    logger_ant.info("Starting read counting...")
    os.makedirs(sample.rcounts)
    task_status, aln_stats = get_readcounts_for_features(assembly_fasta=sample.inputs["assemblyFasta"], 
                                                         gff_file=sample.outputs["annotationGff"], 
                                                         reads_dir=sample.reads,
                                                         bins_dir=sample.bins,
                                                         work_dir=sample.swd, 
                                                         output_dir=sample.rcounts,
                                                         threads=threads,
                                                         bin_extension="fasta",
                                                         feature_type="CDS", 
                                                         feature_id="ID", 
                                                         min_aln_q=20)
    if task_status:
        status = "FAILED" if task_status == 1 else "FINISHED WITH WARNINGS"
        logger_ant.critical(f"Read counting {status}!")
        logger_ant.critical("Terminating sample processing.")
        sample.set_process_status("readCounting", status)
        sample.status = "Annotation failed"
        return sample
    logger_ant.info("Read counting COMPLETED.")
    for k in aln_stats.keys():
        sample.memoize_metric(k, aln_stats[k])
    sample.set_process_status("readCounting", "COMPLETED")
    sample.status = "Annotated"

    if clean_up:
        logger_ant.info("Removing temporary files...")
        tmps = [sample.swd, sample.reads]
        for tmp in tmps:
            utils.remove_tmp(tmp)

    return sample
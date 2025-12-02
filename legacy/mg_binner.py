import os
import subprocess
import concurrent.futures
import threading
import logging
import re
import pandas as pd
import random
import json

import utils

logger_bin = logging.getLogger("default")

def generate_features_semibin2(concat_contigs, aln_dir, output_dir, threads, 
                               min_contig_len=1000):
    '''
    Запускает функцию generate_sequence_features_multi пакета SemiBin2. 
    Производит генерацию фич для последующего обучения модели 
    на основе последовательностей контигов и глубины их покрытия ридами.
    Работает на наборе из нескольких образцов.

    Аргументы:
        concat_contigs (str)  путь к общему FASTA файлу, включающему контиги всех образцов из набора.
        aln_dir (str)         путь к директории с .bam файлами выравниваниями ридов каждого отдельного образца на
                              общий FASTA файл, включающий контиги всех образцов из набора.
        output_dir (str)      путь к директории для сохранения результатов.
        threads (int)         число потоков.
        min_contig_len (int)  минимальная длина контига (от 1000 п.о.) для его включения в процесс биннинга [1000].
    '''
    env = "binning"
    conda_run_cmd = f"conda run -n {env}"
    semibin2 = "SemiBin2"

    cmd = f"""{semibin2} generate_sequence_features_multi 
              --threads {threads} 
              --min-len {min_contig_len} 
              -i {concat_contigs} 
              -b {os.path.join(aln_dir, "*.bam")} 
              -o {output_dir}"""
    cmd = " ".join(cmd.split())
    logger_bin.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_bin.error("An error occurred!", exc_info=True)
        return 1

def training_model_semibin2(output_dir, sample_id, threads, seed=42):
    '''
    Запускает функцию train_self пакета SemiBin2. Обучает модель.
    Работает на наборе из нескольких образцов.

    Аргументы:
        output_dir (str)  путь к директории для сохранения результатов.
        sample_id (str)   идентификатор образца.
        threads (int)     число потоков.
        seed (int)        seed [42].
    '''
    env = "binning"
    conda_run_cmd = f"conda run -n {env}"
    semibin2 = "SemiBin2"

    cmd = f"""{semibin2} train_self 
              --threads {threads} --random-seed {seed} --engine cpu 
              --data {os.path.join(output_dir, "samples", sample_id, "data.csv")} 
              --data-split {os.path.join(output_dir, "samples", sample_id, "data_split.csv")} 
              --output {os.path.join(output_dir, "models", sample_id)}"""
    cmd = " ".join(cmd.split())
    logger_bin.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_bin.error("An error occurred!", exc_info=True)
        return 1

def bin_semibin2(output_dir, sample_id, contigs_fasta, threads, 
                 min_contig_len=1000, min_bin_size=50, seed=42):
    '''
    Запускает функцию bin_short пакета SemiBin2. Осуществляет биннинг контигов.
    Работает на наборе из нескольких образцов.

    Аргументы:
        output_dir (str)      путь к директории для сохранения результатов.
        sample_id (str)       идентификатор образца.
        contigs_fasta (str)   путь к FASTA файлу метагеномной сборки образца.
        threads (int)         число потоков.
        min_contig_len (int)  минимальная длина контига (от 1000 п.о.) для его включения в процесс биннинга [1000].
        min_bin_size (int)    минимальный размер бина (от 50 kb) [50].
        seed (int)            seed [42].
    '''
    env = "binning"
    conda_run_cmd = f"conda run -n {env}"
    semibin2 = "SemiBin2"

    cmd = f"""{semibin2} bin_short 
              --threads {threads} --random-seed {seed} --engine cpu 
              --min-len {min_contig_len} --minfasta-kbs {min_bin_size} 
              -i {contigs_fasta} 
              --model {os.path.join(output_dir, "models", sample_id, "model.pt")} 
              --data {os.path.join(output_dir, "samples", sample_id, "data.csv")} 
              -o {os.path.join(output_dir, f"{sample_id}_bins")}"""
    cmd = " ".join(cmd.split())
    logger_bin.debug(f"CMD: {cmd}")
    try:
        subprocess.run(f"{conda_run_cmd} {cmd}", shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_bin.error("An error occurred!", exc_info=True)
        return 1

def bin_mtr(output_dir, sample_set, n_workers, threads_per_worker,
            min_contig_len=1000, min_bin_size=50, seed=42):
    '''
    Запускает этапы обучения модели (training_model_semibin2()) и биннинга контигов (bin_semibin2()) 
    разных образцов в параллельных потоках выполнения. 
    
    Использует:
        - concurrent.futures
        - threading
    
    Аргументы:
        output_dir (str)          путь к дирректории, куда будут сохраняться результаты работы.
        sample_set (list)         список экземпляров класса BinnedSample.
        n_workers (int)           число потоков выполнения.
        threads_per_worker (int)  число потоков CPU на поток выполнения.
        min_contig_len (int)      минимальная длина контига (от 1000 п.о.) для его включения в 
                                  процесс биннинга [1000].
        min_bin_size (int)        минимальный размер бина (от 50 kb) [50].
        seed (int)                seed [42].
    '''
    def binning_task(sample, cancellation_event):
        logger_bin.info(f"Sample {sample.id}: running...")

        if cancellation_event.is_set():
            logger_bin.info(f"Sample {sample.id}: binning CANCELLED.")
            return
            
        logger_bin.info(f"Sample {sample.id}: training model...")
        task_status = training_model_semibin2(output_dir=output_dir, 
                                              sample_id=sample.id, 
                                              threads=threads_per_worker, 
                                              seed=seed)
        if task_status:
            raise ValueError(f"Sample {sample.id}: an error occurred while training model!")
            logger_bin.info(f"Sample {sample.id}: binning FAILED.")
            return

        if cancellation_event.is_set():
            logger_bin.info(f"Sample {sample.id}: binning CANCELLED.")
            return
        
        logger_bin.info(f"Sample {sample.id}: binning...")
        task_status = bin_semibin2(output_dir=output_dir, 
                                   sample_id=sample.id, 
                                   contigs_fasta=sample.inputs["assemblyFasta"],
                                   threads=threads_per_worker,
                                   min_contig_len=min_contig_len, 
                                   min_bin_size=min_bin_size, 
                                   seed=seed)
        if task_status:
            raise ValueError(f"Sample {sample.id}: an error occurred while binning!")
            logger_bin.info(f"Sample {sample.id}: binning FAILED.")
            return
        logger_bin.info(f"Sample {sample.id}: binning COMPLETE.")
        return

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=n_workers)
    cancellation_event = threading.Event()
    futures = []
    for sample in sample_set:
        futures.append(executor.submit(binning_task, sample, cancellation_event))

    exception_occurred = False
    done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_EXCEPTION)
    for future in done:
        if future.exception():
            logger_bin.error(future.exception())
            exception_occurred = True
            break
    
    if exception_occurred:
        logger_bin.critical("Cancelling pending and running binning tasks...")
        cancellation_event.set()
        executor.shutdown(wait=True, cancel_futures=True)
        logger_bin.info("All running binning tasks have been completed.")
        logger_bin.info("All pending binning tasks have been cancelled.")
        return 1
    else:
        logger_bin.info("All binning tasks completed successfully.")
        executor.shutdown(wait=True)
        return 0

def process_contig2bin_info(contig_bins_tsv, recluster_bins_info_tsv, output_file, sample_id):
    '''
    Принимает файлы contig_bins.tsv и recluster_bins_info.tsv, генерируемые SemiBin2 для каждого образца 
    в процессе биннинга. Записывает .tsv файл с информацией об отношениях контиг - бин. 
    Дополнителньо возвращает Pandas DataFrame с аналогичным содержанием.
    В процессе работы функции имена бинов приводятся к виду "bin_###".
    Контиги, которые не были отнесены к какому-либо бину, помечаются как принадлежащие бину "bin_000"

    Использует:
        - re
        - pandas
        
    Аргументы:
        contig_bins_tsv (str)          путь к файлу contig_bins.tsv.
        recluster_bins_info_tsv (str)  путь к файлу recluster_bins_info.tsv.
        output_file (str)              итоговый .tsv файл.
        sample_id (str)                идентификатор образца.
    '''
    
    for f in [contig_bins_tsv, recluster_bins_info_tsv]:
        if os.path.exists(f):
           logger_bin.info(f"Sample {sample_id}: file '{f}' found.")
        else:
            logger_ant.error(f"Sample {sample_id}: the required '{f}' does not exist.")
            return 1, None
    
    try:
        bins_info = pd.read_table(recluster_bins_info_tsv, sep="\t")
    except Exception as error:
        logger_bin.error(f"Sample {sample_id}: an error occurred!", exc_info=True)
        return 1, None
    contigs_in_final_bins = [f'bin_{re.search("(?<=SemiBin_)[0-9]+", f).group(0)}' for f in bins_info["filename"]]
    
    try:
        contig2bin = pd.read_table(contig_bins_tsv, sep="\t")
    except Exception as error:
        logger_bin.error(f"Sample {sample_id}: an error occurred!", exc_info=True)
        return 1, None
    contig2bin["bin"] = [f'bin_{n}' for n in contig2bin["bin"]]
    contig2bin.loc[~contig2bin["bin"].isin(contigs_in_final_bins), "bin"] = "bin_000"
    
    try:
        contig2bin.to_csv(output_file, sep='\t', index=False)
        logger_bin.info(f"Sample {sample_id}: contig to bin relations were written to a file: '{output_file}'")
    except Exception as error:
        logger_bin.error(f"Sample {sample_id}: an error occurred while writing contig to bin relations to the '{output_file}'", 
                         exc_info=True)
    return 0, contig2bin

def restore_bin(contigs_fasta, bin_id, contig2bin, output_file, sample_id, 
                formate=False):
    '''
    Генерирует FASTA файл бина из контигов и информации об отношениях контиг - бин, 
    переданной в виде Pandas DataFrame.

    Использует:
        - pandas
        - utils.py
        
    Аргументы:
        contigs_fasta (str)  FASTA файл с контигами (сборка метагенома).
        bin_id (str)         идентификатоб бина вида "bin_###".
        contig2bin           Pandas DataFrame с информацией об отношениях контиг - бин.
        output_file (str)    итоговый FASTA файл.
        sample_id (str)      идентификатор образца.
        formate (bool)       привести FASTA ID к виду "bin_id:contig_id"? [False]
    '''

    contigs = f"{output_file}.contigs"
    tmp_file = f"{output_file}.tmp"
    try:
        with open(contigs, "x") as o_file:
            o_file.write("\n".join(list(contig2bin.loc[contig2bin["bin"] == bin_id, "contig"])))
    except Exception as error:
        logger_bin.error("An error occurred!", exc_info=True)
        return 1

    task_status = utils.run_seqkit_grep(input_file=contigs_fasta, 
                                        patterns_file=contigs,
                                        output_file=tmp_file if formate else output_file)
    if task_status:
        logger_bin.warning(f"Sample {sample_id}: unable to restore bin {bin_id} from assembly .fasta file.")
        return 1

    if formate:
        task_status = utils.run_seqkit_replace(input_file=tmp_file, 
                                               output_file=output_file, 
                                               pattern="^", 
                                               replacement=f"{bin_id}:")
        if task_status:
            logger_bin.warning(f"Sample {sample_id}: unable to formate {bin_id} .fasta file.")
            return 1
    try:
        os.remove(contigs)
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
    except Exception as error:
        logger_bin.error("An error occurred!", exc_info=True)
    return 0

def replace_assembly_fasta(assembly_fasta, bins_dir, ext="fasta"):
    '''
    Заменяет оригинальный FASTA файл со сборкой метагенома новым, скомбинированным из
    FASTA файлов бинов. Важно, поскольку в FASTA файлах бинов заголовки записей уже приведены к 
    виду "bin_###:contig_###". Далее в процессе аннотации используются идентификаторы контигов
    именно такого формата.
    ______________________
    Заменять файл тупо, но, как оказалось, более эффективно по времени, 
    нежели редактирование каждой записи.
        
    Аргументы:
        assembly_fasta (str)  путь к оригинальному FASTA файл с контигами (сборка метагенома).
        bins_dir (str)        путь к директории с FASTA файлами бинов.
        ext (str)             расширение FASTA файлов ["fasta"].
    '''
    try:
        os.remove(assembly_fasta)
    except Exception as error:
        logger_bin.error("An error occurred!", exc_info=True)
        logger_bin.warning("Unable to remove original assembly .fasta file.")
        return 1

    cmd = f"cat {bins_dir}/*.{ext} > {assembly_fasta}"
    logger_bin.debug(f"CMD: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        return 0
    except subprocess.CalledProcessError as error:
        logger_bin.error("An error occurred!", exc_info=True)
        return 1

def bins_post_processing_mtr(binning_dir, sample_set, n_workers):
    '''
    Запускает постпроцессинг результатов бининга разных образцов в параллельных потоках выполнения. 
    
    Использует:
        - concurrent.futures
        - threading
    
    Аргументы:
        binning_dir (str)         путь к директории с результатами биннинга набора образцов.
        sample_set (list)         список экземпляров класса BinnedSample.
        n_workers (int)           число потоков выполнения.
    '''
    def bins_pp_task(sample):
        logger_bin.info(f"Sample {sample.id}: processing bins info...")
        recluster_bins_info_tsv = os.path.join(binning_dir, f"{sample.id}_bins", "recluster_bins_info.tsv")
        contig_bins_tsv = os.path.join(binning_dir, f"{sample.id}_bins", "contig_bins.tsv")
        task_status, contig2bin = process_contig2bin_info(contig_bins_tsv=contig_bins_tsv, 
                                                          recluster_bins_info_tsv=recluster_bins_info_tsv, 
                                                          output_file=sample.outputs["contigToBin"], 
                                                          sample_id=sample.id)
        if task_status:
            logger_bin.warning(f"Sample {sample.id}: an error occurred while processing bins info.")
            return 1, sample.id
        
        logger_bin.info(f"Sample {sample.id}: processing bins .fasta files...")
        got_warnings = False

        os.makedirs(sample.bins)
        semibin2_bins_dir = os.path.join(binning_dir, f"{sample.id}_bins", "output_bins")
        mg_bins = os.listdir(semibin2_bins_dir)
        for mg_bin in mg_bins:
            mg_bin_new = f"bin_{re.search("(?<=SemiBin_)[0-9]+", mg_bin).group(0)}"
            task_status = utils.run_seqkit_replace(input_file=os.path.join(semibin2_bins_dir, mg_bin), 
                                                   output_file=os.path.join(sample.bins, f"{mg_bin_new}.fasta"), 
                                                   pattern="^", 
                                                   replacement=f"{mg_bin_new}:")
            if task_status:
                logger_bin.warning(f"Sample {sample.id}: an error occurred while processing {mg_bin}.")
                return 1, sample.id
                
        logger_bin.info(f"Sample {sample.id}: merging unbinned contigs into bin_000...")
        task_status = restore_bin(contigs_fasta=sample.inputs["assemblyFasta"], 
                                  bin_id="bin_000", 
                                  contig2bin=contig2bin, 
                                  output_file=os.path.join(sample.bins, "bin_000.fasta"), 
                                  sample_id=sample.id,
                                  formate=True)
        if task_status:
            logger_bin.warning(f"Sample {sample.id}: merging unbinned contigs into bin_000 FAILED!")
            return 1, sample.id

        logger_bin.info(f"Sample {sample.id}: replacing original assembly .fasta file...")
        task_status = replace_assembly_fasta(assembly_fasta=sample.inputs["assemblyFasta"], 
                                             bins_dir=sample.bins, 
                                             ext="fasta")
        if task_status:
            logger_bin.warning(f"Sample {sample.id}: replacing original assembly .fasta file FAILED!")
            return 2, sample.id
        
        logger_bin.info(f"Sample {sample.id}: {len(mg_bins) + 1} bins .fasta files have been saved to '{sample.bins}'.")
        return 0, sample.id

    logger_bin.debug(f"Runnig in {n_workers} threads.")
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for sample in sample_set:
            futures.append(executor.submit(bins_pp_task, sample))

        done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.ALL_COMPLETED)
        
        save_semibin2_outputs = False
        for future in done:
            try:
                result = future.result()
                if result[0] == 1:
                    logger_bin.warning(f"Future {future}: sample {result[1]} terminated!")
                    save_semibin2_outputs = True
                elif result[0] == 2:
                    logger_bin.warning(f"Future {future}: sample {result[1]} finished with warnings!")
                    save_semibin2_outputs = True
                else:
                    logger_bin.info(f"Future {future}: sample {result[1]} finished successfully.")
            except Exception as error:
                logger_bin.error(f"An error occurred!", exc_info=True)
                save_semibin2_outputs = True

    if save_semibin2_outputs:
        return 1
    return 0

def bin_mg(sample_set, binning_dir, threads, 
           threads_per_worker=4, min_contig_len=1000, 
           min_bin_size=50, seed=42, clean_up=True):
    '''
    Реализует полный процесс биннига метагеномных сборок с использованием пайплайна Semibin2: 
    0)  Проверяет существование требуемых файлов.
    1)  Объедиянет контиги метагеномных сборок в общий FASTA файл.
    2)  Индексирует полученный объединенный FASTA файл и производит выравнивание на него ридов 
        каждого образца по отдельности. Выравниваются парные и потерявшие пару риды, прошедшие тримминг
        (перед сборкой). Выполняет постпроцессинг выравниваний с использованием функционала samtools.
    3)  Запускает пайплайн Semibin2. Этапы обучения модели и биннинга выполняются в несколько програмных потоков.
    4)  Выполняет построцессинг результатов биннинга: генерирует файл 'contig2bin.tsv' с информацией о
        отношениях контиг - бин; объединяет в bin_000 контиги, не вошедшие в финальные бины, 
        и генерирует bin_000.fasta; переформатирует FASTA файлы бинов и метагеномных сборок образцов;
        FASTA файлы бинов переносятся в директорию ./assembly/bins.
    В процессе работы собирается информация об общем числе бинов, полученных для конкретного образца.

    В результате работы в директории сборки образца остаются:
    sample_dir
    |----- assembly.log
    |----- assembly
    |      |----- contigs.fasta
    |      |----- assembly_graph.fastg
    |      |----- contigs.paths
    |      |----- scaffolds.paths
    |      |----- contig2bin.tsv
    |      |----- bins
    |             |----- bin_000.fasta
    |             |----- bin_0.fasta
    |             |----- ...
    |----- reports
           |----- spades.log
           |----- fastp_report.html
           |----- fastp_report.json
    
    Использует:
        - sample_cl.py
        - utils.py
        
    Аргументы:
        sample_set (list)         список экземпляров класса BinnedSample, определяющий структуру содержания
                                  директории сборки образца и хранящий информацию о статусе биннига, 
                                  файлах с важными результатами, а также собираемые метрики.
        binning_dir (str)         путь к директории для сохранения промежуточных результатов биннинга.
                                  Не будет удалена, если в процессе биннига возникнуть некритические 
                                  ошибки, либо clean_up=False.
        threads (int)             число потоков CPU.
        threads_per_worker (int)  число потоков CPU на поток выполнения. 
        min_contig_len (int)      минимальная длина контига (от 1000 п.о.) для его включения в 
                                  процесс биннинга [1000].
        min_bin_size (int)        минимальный размер бина (от 50 kb) [50].
        seed (int)                seed [42].
        clean_up (bool)           удалить все промежуточные файлы? [True] 

    Возвращает:
        - статус сборки образца ["Binned" | "Binned (warnings)" | "Binning failed"]
    '''
    logger_bin.debug(f"Working in '{os.getcwd()}'.")
    
    logger_bin.info("Verifying inputs...")
    all_passed = True
    for sample in sample_set:
        try:
            sample.check_inputs()
            sample.set_process_status("inputsCheck", "COMPLETED")
        except FileNotFoundError as error:
            logger_bin.error(f"Sample {sample.id}: inputs verification FAILED!", exc_info=True)
            sample.set_process_status("inputsCheck", "FAILED")
            sample.status = "Binning failed"
            all_passed = False
    if all_passed:
        logger_bin.info("Inputs verification COMPLETED.")
    else:
        logger_bin.critical("Some samples did not pass the inputs check.")
        logger_bin.critical("Terminating binning.")
        return sample_set
    
    logger_bin.info("Concatenating contigs...")
    concat_contigs_fasta = os.path.join(binning_dir, "concat_contigs.fasta")
    all_passed = True
    c = 1
    for sample in sample_set:
        task_status = utils.run_seqkit_replace(input_file=sample.inputs["assemblyFasta"], 
                                               output_file=concat_contigs_fasta, 
                                               pattern="^", 
                                               replacement=f"{sample.id}:")
        if task_status:
            logger_bin.critical(f"Addition {sample.id} sample contigs FAILED!")
            sample.set_process_status("aligning", "FAILED")
            sample.status = "Binning failed"
            all_passed = False
        else:
            logger_bin.info(f"[{c}] Sample {sample.id} was added.")
        c += 1
    if all_passed:
        logger_bin.info("Contigs concatenation COMPLETED.")
    else:
        logger_bin.critical("Some samples contigs were not concateneted.")
        logger_bin.critical("Terminating binning.")
        return sample_set

    logger_bin.info("Aligning samples reads on concatenated contigs...")
    aln_dir = os.path.join(binning_dir, "alignments")
    index_dir = os.path.join(aln_dir, "index")
    os.makedirs(index_dir)
    logger_bin.info("Building index...")
    task_status = utils.run_bowtie2_build(fasta_file=concat_contigs_fasta, 
                                          output_dir=index_dir, 
                                          threads=threads)
    if task_status:
        logger_bin.critical("Building index FAILED!")
        logger_bin.critical("Terminating binning.")
        for sample in sample_set:
            sample.set_process_status("aligning", "FAILED")
            sample.status = "Binning failed"
        return sample_set
    logger_bin.info("Building index COMPLETED.")

    logger_bin.info("Performing per sample alignment...")
    all_passed = True
    c = 1
    for sample in sample_set:
        logger_bin.info(f"[{c}] Sample {sample.id}: aligning...")
        
        sam_file=os.path.join(aln_dir, f"{sample.id}.sam")
        unsorted_bam_file = os.path.join(aln_dir, f"{sample.id}_unsorted.bam")
        bam_file = os.path.join(aln_dir, f"{sample.id}.bam")
        
        task_status = utils.run_bowtie2(reads_dir=sample.reads, 
                                        index_path=os.path.join(index_dir, "index"), 
                                        output_file=os.path.join(aln_dir, f"{sample.id}.sam"), 
                                        threads=threads,
                                        use_unpaired=True)
        if task_status:
            logger_bin.critical(f"[{c}] Sample {sample.id}: aligning FAILED!")
            logger_bin.critical("Terminating binning.")
            sample.set_process_status("aligning", "FAILED")
            sample.status = "Binning failed"
            return sample_set
        logger_bin.info(f"[{c}] Sample {sample.id}: converting .sam to .bam and filtering unmapped reads...")
        task_status = utils.convert_sam_to_bam(sam_file=sam_file, 
                                               bam_file=unsorted_bam_file, 
                                               threads=threads,
                                               exclude_unmapped=True)
        if task_status:
            logger_bin.critical(f"[{c}] Sample {sample.id}: converting .sam to .bam and filtering unmapped reads FAILED!")
            logger_bin.critical("Terminating binning.")
            sample.set_process_status("aligning", "FAILED")
            sample.status = "Binning failed"
            return sample_set
        logger_bin.info(f"[{c}] Sample {sample.id}: sorting .bam...")
        task_status = utils.sort_bam(unsorted_bam_file=unsorted_bam_file, 
                                     sorted_bam_file=bam_file, 
                                     threads=threads)
        if task_status:
            logger_bin.critical(f"[{c}] Sample {sample.id}: sorting .bam FAILED!")
            logger_bin.critical("Terminating binning.")
            sample.set_process_status("aligning", "FAILED")
            sample.status = "Binning failed"
            return sample_set
        logger_bin.info(f"[{c}] Sample {sample.id}: removing .sam and unsorted .bam files...")
        try:
            os.remove(sam_file)
            os.remove(unsorted_bam_file)
        except exception as error:
            logger_bin.error(f"[{c}] Sample {sample.id}: an error occurred!", exc_info=True)
            logger_bin.critical(f"[{c}] Sample {sample.id}: unable to remove .sam and unsorted .bam files.")
            logger_bin.critical(f"[{c}] Sample {sample.id}: this will affect binning process!")
            logger_bin.critical("Terminating binning.")
            sample.set_process_status("aligning", "FAILED")
            sample.status = "Binning failed"
            return sample_set
        sample.set_process_status("aligning", "COMPLETED")
        logger_bin.info(f"[{c}] Sample {sample.id}: COMPLETED.")
        c += 1
    logger_bin.info("Per sample aligning COMPLETED.")

    logger_bin.info("Binning samples with SemiBin2...")
    logger_bin.info("Generation of sequence features...")
    task_status = generate_features_semibin2(concat_contigs=concat_contigs_fasta, 
                                             aln_dir=aln_dir, 
                                             output_dir=binning_dir, 
                                             threads=threads, 
                                             min_contig_len=min_contig_len)
    if task_status:
        logger_bin.critical("Generation of sequence features FAILED!")
        logger_bin.critical("Terminating binning.")
        for sample in sample_set:
            sample.set_process_status("binning", "FAILED")
            sample.status = "Binning failed"
        return sample_set
    logger_bin.info("Generation of sequence features COMPLETED.")
    
    logger_bin.info("Performing per sample model training and binning...")
    #### MOCK_POOL ####
    real_samples = [sample for sample in sample_set if sample.pool != "MOCK_POOL"]
    sample_set = real_samples
    ###################
    n_tasks = len(sample_set)
    threads_per_task = threads // n_tasks
    if threads_per_task > threads_per_worker:
        threads_per_worker = threads_per_task
    max_workers = threads // threads_per_worker
    max_workers = max_workers if max_workers < n_tasks else n_tasks
    logger_bin.debug(f"Workers: {max_workers}")
    logger_bin.debug(f"Threads per worker: {threads_per_worker}")
    task_status = bin_mtr(output_dir=binning_dir, 
                          sample_set=sample_set, 
                          min_contig_len=min_contig_len, 
                          min_bin_size=min_bin_size, 
                          n_workers=max_workers, 
                          threads_per_worker=threads_per_worker, 
                          seed=seed)
    if task_status:
        logger_bin.critical("Performing per sample model training and binning FAILED!")
        logger_bin.critical("Terminating binning.")
        for sample in sample_set:
            sample.set_process_status("binning", "FAILED")
            sample.status = "Binning failed"
        return sample_set
    for sample in sample_set:
            sample.set_process_status("binning", "COMPLETED")
    logger_bin.info("Performing per sample model training and binning COMPLETED.")

    logger_bin.info("Bins post-processing...")
    n_tasks = len(sample_set)
    max_workers = threads if threads < n_tasks else n_tasks
    logger_bin.debug(f"Workers: {max_workers}")
    task_status = bins_post_processing_mtr(binning_dir=binning_dir, 
                                           sample_set=sample_set, 
                                           n_workers=max_workers)
    if task_status:
        logger_bin.warning("Bins post-processing finished with warnings.")
        logger_bin.info(f"The intermediate binning files will be saved and can be found in the '{binning_dir}'.")
        for sample in sample_set:
            sample.set_process_status("binsProcessing", "FAILED")
            sample.status = "Binned (warnings)"
        return sample_set
    for sample in sample_set:
        sample.set_process_status("binsProcessing", "COMPLETED")
        sample.status = "Binned"
    logger_bin.info("Bins post-processing COMPLETED.")

    logger_bin.info("Сollecting information about the number of bins per sample...")
    c = 1
    for sample in sample_set:
        contig2bin = pd.read_csv(sample.outputs["contigToBin"], sep="\t")
        n_bins = contig2bin["bin"].nunique()
        sample.memoize_metric("nBins", n_bins)
        logger_bin.info(f"[{c}] Sample {sample.id}: {n_bins} bins.")
        c += 1
        
    if clean_up:
        logger_bin.info("Removing temporary files...")
        utils.remove_tmp(binning_dir)
            
    return sample_set
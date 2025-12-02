import os
import argparse
import logging
import logging.config
import json
import random
import string
import time
from datetime import datetime

from sample_cl import AssembledSample, BinnedSample, AnnotatedSample
from mg_assembler import assemble_mg
from mg_binner import bin_mg
from mg_annotator import annotate_mg

def get_logger_config(log_files, mode="w", level="DEBUG"):
    config = {
        "disable_existing_loggers": False,
        "version" : 1,
        "formatters" : {
            "default" : {
                "format" : "%(asctime)s [%(filename)s:%(funcName)s] %(levelname)s: %(message)s",
                "datefmt" : "%Y-%m-%d %H:%M:%S",
            }
        },
        "handlers" : {
            "stdout" : {
                "class" : "logging.StreamHandler",
                "level" : level,
                "formatter" : "default",
                "stream" : "ext://sys.stdout"
            }
        },
        "loggers" : {
            "default" : {
                "handlers" : [],
                "level" : level,
                "propagate" : False
            }
        }
    }

    for i in range(len(log_files)):
        config["handlers"][f"log_file_{i}"] = {
            "class" : "logging.FileHandler",
            "level" : level,
            "formatter" : "default",
            "filename" : log_files[i],
            "mode" : mode
        }
    config["loggers"]["default"]["handlers"] = list(config["handlers"].keys())
    return config

def add_header_to_log(log_file, header):
    header = "\n  " + "=" * (len(header) + 4) + "\n" + f"    {header.upper()}" + "\n  " + "="  * (len(header) + 4) + "\n\n"
    with open(log_file, "a") as log:
        log.write(header)

def generate_rstr(str_len):
    characters = string.ascii_letters + string.digits 
    rstr = ''.join(random.choices(characters, k=str_len))
    return rstr

def formate_et(elapsed_sec):
    h = int(elapsed_sec // 3600)
    m = int((elapsed_sec % 3600) // 60)
    s = int(elapsed_sec % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

def run_assembly(args):
    CURRENT_DIR = os.getcwd()
    os.chdir(args.output_dir)
    START = time.time()
    
    pool_id, sample_id = args.mgid.split(":")
    fastqs = {"R1" : [os.path.join(args.fastq_dir, fq) for fq in args.r1_fastqs],
              "R2" : [os.path.join(args.fastq_dir, fq) for fq in args.r2_fastqs]}

    sample = AssembledSample(sample_id, pool_id, fastqs)
    os.makedirs(sample.dir, exist_ok=True)
    logging.config.dictConfig(get_logger_config(log_files=[sample.log_file], mode="a"))
    add_header_to_log(sample.log_file, "assembly")
        
    sample = assemble_mg(sample=sample, 
                         threads=args.threads, 
                         memory_limit=args.memory_limit, 
                         min_q=args.min_q, 
                         read_min_len=args.read_min_len, 
                         deduplication=args.disable_dedup,
                         contig_min_len=args.contig_min_len, 
                         contig_basename=args.contig_basename, 
                         clean_up=args.disable_cleaning)

    _sample_info = sample.formate()
    sample_info = {
        "mgid" : _sample_info["mgid"], 
        "sampleStatus" : _sample_info["sampleStatus"],
        "processStatus" : {
            "assembly" : _sample_info["processStatus"]
        },
        "outputs" : _sample_info["outputs"],
        "metrics" : _sample_info["metrics"]
    }
    with open(os.path.join(sample.dir, "sample.json"), "w") as s_info_json:
        json.dump(sample_info, s_info_json, indent=4)

    wall_time = formate_et(time.time() - START)
    dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open("sample_processing.log", "a") as sp_log:
        sp_log.write(f"{dt} | Sample {args.mgid}: {sample.status} [{wall_time}].\n")
    os.chdir(CURRENT_DIR)
    return

def run_binning(args):
    CURRENT_DIR = os.getcwd()
    os.chdir(args.output_dir)

    START = time.time()
    ss_id = f"SS_{generate_rstr(9)}"
    if os.path.exists("sample_sets.json"):
        with open("sample_sets.json", "r") as ss_json:
            ss = json.load(ss_json)
        ss[ss_id] = args.mgids
    else:
        ss = {ss_id : args.mgids}
    with open("sample_sets.json", "w") as ss_json:
        json.dump(ss, ss_json, indent=4)
        
    sample_set = []
    for mgid in args.mgids:
        pool_id, sample_id = mgid.split(":")
        sample_set.append(BinnedSample(sample_id, pool_id))

    log_files = [sample.log_file for sample in sample_set if sample.pool != "MOCK_POOL"] #### MOCK_POOL ####
    logging.config.dictConfig(get_logger_config(log_files=log_files, mode="a"))
    for f in log_files:
        add_header_to_log(f, "binning")
        
    os.makedirs(ss_id, exist_ok=False)
    sample_set = bin_mg(sample_set=sample_set, 
                        binning_dir=ss_id, 
                        threads=args.threads, 
                        threads_per_worker=args.threads_per_worker, 
                        min_contig_len=args.contig_min_len, 
                        min_bin_size=args.min_bin_size, 
                        seed=42, 
                        clean_up=args.disable_cleaning)
    
    wall_time = formate_et(time.time() - START)
    dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for sample in sample_set:
        if sample.pool != "MOCK_POOL":  #### MOCK_POOL ####
            _sample_info = sample.formate()
            with open(os.path.join(sample.dir, "sample.json"), "r") as s_info_json:
                json_info = json.load(s_info_json)
            sample_info = {
                "mgid" : _sample_info["mgid"], 
                "sampleStatus" : _sample_info["sampleStatus"],
                "processStatus" : {
                    **json_info["processStatus"],
                    "binning" : _sample_info["processStatus"]
                },
                "outputs" : {
                    **json_info["outputs"],
                    **_sample_info["outputs"],
                },
                "metrics" : {
                    **json_info["metrics"],
                    **_sample_info["metrics"],
                }
            }
            with open(os.path.join(sample.dir, "sample.json"), "w") as s_info_json:
                json.dump(sample_info, s_info_json, indent=4)
        with open("sample_processing.log", "a") as sp_log:
            sp_log.write(f"{dt} | Sample {sample.pool}:{sample.id}: {sample.status} [{wall_time}].\n")

    os.chdir(CURRENT_DIR)
    return

def run_annotation(args):
    os.environ["GTDBTK_DATA_PATH"] = args.gtdb
    os.environ["CHECKM_DATA_PATH"] = args.checkm_db
    os.environ["CHECKM2DB"] = args.checkm2_db
    
    CURRENT_DIR = os.getcwd()
    os.chdir(args.output_dir)

    START = time.time()
    
    pool_id, sample_id = args.mgid.split(":")
    sample = AnnotatedSample(sample_id, pool_id)
    
    logging.config.dictConfig(get_logger_config(log_files=[sample.log_file], mode="a"))
    add_header_to_log(sample.log_file, "annotation")

    sample = annotate_mg(sample=sample,
                         ref_db=args.mmseqs2_db,
                         threads=args.threads, 
                         use_aa_seqs_for_bins_qc=args.disable_qc_with_aa, 
                         checkm_v_for_mags_q_asses=args.checkm_v,
                         save_cds_nuc_seqs=args.dont_save_cds_nuc_seqs,
                         use_repr_seqs_for_cds_annot=args.disable_cds_clust,
                         clean_up=args.disable_cleaning)

    _sample_info = sample.formate()
    with open(os.path.join(sample.dir, "sample.json"), "r") as s_info_json:
        json_info = json.load(s_info_json)
    sample_info = {
        "mgid" : _sample_info["mgid"], 
        "sampleStatus" : _sample_info["sampleStatus"],
        "processStatus" : {
            **json_info["processStatus"],
            "annotation" : _sample_info["processStatus"]
        },
        "outputs" : {
            **json_info["outputs"],
            **_sample_info["outputs"],
        },
        "metrics" : {
            **json_info["metrics"],
            **_sample_info["metrics"],
        }
    }
    with open(os.path.join(sample.dir, "sample.json"), "w") as s_info_json:
        json.dump(sample_info, s_info_json, indent=4)

    wall_time = formate_et(time.time() - START)
    dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open("sample_processing.log", "a") as sp_log:
        sp_log.write(f"{dt} | Sample {args.mgid}: {sample.status} [{wall_time}].\n")
    os.chdir(CURRENT_DIR)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ассемблер метагеномных шотганов.")
    subparsers = parser.add_subparsers(dest="command", help="Доступный функционал:")

    run_asmb_parser = subparsers.add_parser("assemble", help="Запускает сборку образца.")
    run_asmb_parser.add_argument("mgid", type=str, help="Идентификатор образца в формате 'pool_id:sample_id'.")
    run_asmb_parser.add_argument("--r1_fastqs", type=str, required=True, nargs="+", help="Список R1 .fastq файлов.")
    run_asmb_parser.add_argument("--r2_fastqs", type=str, required=True, nargs="+", help="Список R2 .fastq файлов.")
    run_asmb_parser.add_argument("-f", "--fastq_dir", type=str, default="/inputs", help="Путь к директории с R1 и R2 .fastq файлами. Если указан, то --r#_fastqs должны передаваться ТОЛЬКО имена файлов. При запуске в Docker использовать '/inputs'.")
    run_asmb_parser.add_argument("-o", "--output_dir", type=str, default="/outputs", help="Путь к директории для хранения резултатов. В ней будет создана директория образца './output_dir/pool_id/sample_id'. При запуске в Docker использовать '/outputs'.")
    run_asmb_parser.add_argument("-T", "--threads", type=int, default=36, help="Число потоков CPU.")
    run_asmb_parser.add_argument("-M", "--memory_limit", type=int, default=300, help="Лимит память в Gb.")
    run_asmb_parser.add_argument("-q", "--min_q", type=int, default=20, help="Минимальное качество основания.")
    run_asmb_parser.add_argument("-l", "--read_min_len", type=int, default=70, help="Минимальная допустимая длина рида после тримминга.")
    run_asmb_parser.add_argument("--disable_dedup", action="store_false", help="Не выполнять дедупликацию?")
    run_asmb_parser.add_argument("-L", "--contig_min_len", type=int, default=1000, help="Минимальная длина контига.")
    run_asmb_parser.add_argument("-c", "--contig_basename", type=str, default="contig", help="Базовое имя контига для использования в FASTA IDs.")
    run_asmb_parser.add_argument("--disable_cleaning", action="store_false", help="Не удалять промежуточные файлы?")
    run_asmb_parser.set_defaults(func=run_assembly)

    run_bin_parser = subparsers.add_parser("bin", help="Запускает биннинг набора образцов.")
    run_bin_parser.add_argument("mgids", type=str, nargs="+", help="Идентификаторы образцов формата 'pool_id:sample_id'.")
    run_bin_parser.add_argument("-o", "--output_dir", type=str, default="/outputs", help="Путь к директории для сохранения резултатов. В ней будет создана директория для временного хранения результатов биннинга './output_dir/ss_id' (ss_id имеет формат 'SS_#########'), которая не будет удалена в случае возникновения некритических ошибок в процессе биннинга или использования флага '--disable_cleaning'. При запуске в Docker использовать '/outputs'.")
    run_bin_parser.add_argument("-T", "--threads", type=int, default=36, help="Число потоков CPU.")
    run_bin_parser.add_argument("-t", "--threads_per_worker", type=int, default=4, help="Минимальное число потоков CPU на поток выполнения (для этапов, реализованных через threading).")
    run_bin_parser.add_argument("-L", "--contig_min_len", type=int, default=1000, help="Минимальная длина контига.")
    run_bin_parser.add_argument("-s", "--min_bin_size", type=int, default=50, help="Минимальный размер бина (в Kbase).")
    run_bin_parser.add_argument("--disable_cleaning", action="store_false", help="Не удалять промежуточные файлы?")
    run_bin_parser.set_defaults(func=run_binning)

    run_ant_parser = subparsers.add_parser("annotate", help="Запускает аннотацию образца.")
    run_ant_parser.add_argument("mgid", type=str, help="Идентификатор образца формата 'pool_id:sample_id'.")
    run_ant_parser.add_argument("--gtdb", type=str, default="/dbs/gtdb/release226", help="Путь к GTDB. При запуске в Docker желательно, чтобы эта ДБ локализовалась в директрии, монтируемой в '/dbs'.")
    run_ant_parser.add_argument("--checkm_db", type=str, default="/dbs/checkmdb", help="Путь к БД CheckM. При запуске в Docker желательно, чтобы эта ДБ локализовалась в директрии, монтируемой в '/dbs'.")
    run_ant_parser.add_argument("--checkm2_db", type=str, default="/dbs/checkm2db/CheckM2_database/uniref100.KO.1.dmnd", help="Путь к БД CheckM2. При запуске в Docker желательно, чтобы эта ДБ локализовалась в директрии, монтируемой в '/dbs'.")
    run_ant_parser.add_argument("--mmseqs2_db", type=str, default="/dbs/mmseqs2dbs/UniRef90", help="Путь к референсной БД для аннотации CDS. Принимаются референсные БД, cобранные под MMseqs2. По умолчанию: UniRef90. При запуске в Docker желательно, чтобы эта ДБ локализовалась в директрии, монтируемой в '/dbs'.")
    run_ant_parser.add_argument("-o", "--output_dir", type=str, default="/outputs", help="Путь к директории для сохранения резултатов. Результаты аннотации будут сохраняться в уже существующую директорию сборки образца './output_dir/pool_id/sample_id'. При запуске в Docker использовать '/outputs'.")
    run_ant_parser.add_argument("-T", "--threads", type=int, default=36, help="Число потоков CPU.")
    run_ant_parser.add_argument("--disable_qc_with_aa", action="store_false", help="Не использовать аминокислотные последовательности CDS для MAGs QC? Будут использованы полные нуклеотидные последовательности.")
    run_ant_parser.add_argument("--checkm_v", type=int, default=1, choices=[1, 2], help="Версия CheckM, результаты QC которой будут использованы для присвоения MAGs категорий качества.")
    run_ant_parser.add_argument("--dont_save_cds_nuc_seqs", action="store_false", help="Не сохранять нуклеотидные последовательности CDS?")
    run_ant_parser.add_argument("--disable_cds_clust", action="store_false", help="Не выполянть кластеризацию последовательности CDS перед аннотацией?")
    run_ant_parser.add_argument("--disable_cleaning", action="store_false", help="Не удалять промежуточные файлы?")
    run_ant_parser.set_defaults(func=run_annotation)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        parser.print_help()

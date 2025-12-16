#!/usr/bin/env python3

import sys
from Bio import SeqIO # type: ignore
from pathlib import Path

def process_contigs(
                    unfiltered_contigs: Path,
                    min_len: int,
                    contig_basename: str,
                    description: bool,
                    additional_meta: str
                   ) -> int:
    """
    Осуществляет постобработку контигов: фильтрацию по длине, переименование и формирование описания.

    Функция читает исходный FASTA-файл, полученный от MetaSPAdes, и записывает в выходной файл
    только те контиги, длина которых не меньше заданного порога. Каждый контиг переименовывается
    с использованием указанного префикса. В описание (description) каждого контига добавляется
    его длина, покрытие k-мерами (извлекается из исходного ID) и опциональная дополнительная информация.

    Использует модуль SeqIO из BioPython для работы с FASTA-файлами.

    :param unfiltered_contigs: Путь к исходному FASTA-файлу, сгенерированному MetaSPAdes.
    :param filtered_contigs: Путь к выходному FASTA-файлу, в который будут записаны отфильтрованные контиги.
    :param min_len: Минимальная длина контига для включения в выходной файл. По умолчанию 1000.
    :param contig_basename: Базовое имя (префикс) для переименования контигов. По умолчанию "contig".
    :param description: Если True, описание (description) контигов будет добавлено. По умолчанию True.
    :param additional_meta: Дополнительная строка, которая будет добавлена к описанию каждого контига. По умолчанию пустая строка.
    :return: Код возврата: 0 при успешном выполнении, 1 при возникновении исключения.
    """
    def get_new_filename(
                         old_filename: Path,
                         old_substring: str,
                         new_substring: str
                        ) -> Path:
        """
        Заменяет подстроки в базовом имени файла
        """
        new_file_basename = old_filename.name.replace(old_substring, new_substring)
        return old_filename.parent / new_file_basename
    
    exit_code = 0
    if unfiltered_contigs.exists(follow_symlinks=True):
        filtered_contigs = get_new_filename(unfiltered_contigs, '.fa', 'processed_contigs.fa')
        try:
            with open(unfiltered_contigs, 'r') as i_file, open(filtered_contigs, 'a') as o_file:
                for contig in SeqIO.parse(i_file, "fasta"):
                    # Фильтрация по длине
                    if len(contig.seq) < min_len:
                        continue

                    # Извлечение информации из исходного ID (формат SPAdes: NODE_1_length_2097_cov_3.0)
                    contig_name_parts = contig.id.split("_")
                    if len(contig_name_parts) < 6:
                        print(f"Некорректный формат ID контига: {contig.id}. Пропускаем.")
                        continue

                    # Переименовываем контиг
                    contig.id = f"{contig_basename}_{contig_name_parts[1]}"

                    # Формируем описание
                    if description:
                        
                        length = len(contig.seq)
                        spades_cov = contig_name_parts[5]
                        contig.description = f"length={length};spades_cov={spades_cov};{additional_meta}".rstrip(';')
                    else:
                        contig.description = ""

                    # Записываем отфильтрованный контиг
                    SeqIO.write(contig, o_file, "fasta")

            return 0

        except FileNotFoundError as e:
            print(f"Файл не найден: {e}")
            return 1
        except PermissionError as e:
            print(f"Ошибка доступа к файлу: {e}")
            return 1
        except Exception as e:
            print(f"Неожиданная ошибка при обработке контигов: {e}")
            return 1
    else:
        print("В рабочей папке не найдены FASTA-файлы")
        exit_code = 127
    return exit_code


if __name__ == "__main__":
    fasta_extensions = ['.fa', '.fasta']

    min_len = int(sys.argv[1])
    contig_basename = sys.argv[2]
    description = bool(sys.argv[3])
    additional_meta = sys.argv[4]
    
    raw_fastas = []
    for ext in fasta_extensions:
        ext_fastas = [f for f in Path('./').rglob(pattern=ext)]
        raw_fastas.extend(ext_fastas)

    for raw_fasta in raw_fastas:
        process_contigs(
                        unfiltered_contigs=raw_fasta,
                        min_len=min_len,
                        contig_basename=contig_basename,
                        description=description,
                        additional_meta=additional_meta
                       )

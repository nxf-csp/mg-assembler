import sys
from Bio import SeqIO
from pathlib import Path

def process_contigs(
    filtered_contigs: str,
    min_len: int = 1000,
    contig_basename: str = "contig",
    no_description: bool = False,
    add_to_description: str = ""
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
    :param no_description: Если True, описание (description) контигов не будет добавлено. По умолчанию False.
    :param add_to_description: Дополнительная строка, которая будет добавлена к описанию каждого контига. По умолчанию пустая строка.
    :return: Код возврата: 0 при успешном выполнении, 1 при возникновении исключения.
    """
    exit_code = 0
    unfiltered_contigs = next(Path('./').resolve().rglob('*.fasta'), None)
    if unfiltered_contigs:
        filtered_contigs = unfiltered_contigs.replace('.fasta', '_filtered.fasta')
        try:
            with open(unfiltered_contigs, 'r') as i_file, open(filtered_contigs, 'a') as o_file:
                for contig in SeqIO.parse(i_file, "fasta"):
                    # Фильтрация по длине
                    if len(contig.seq) < min_len:
                        continue

                    # Извлечение информации из исходного ID (формат SPAdes: NODE_1_length_2097_cov_3.0)
                    contig_name_parts = contig.id.split("_")
                    if len(contig_name_parts) < 6:
                        logger_asmb.warning(f"Некорректный формат ID контига: {contig.id}. Пропускаем.")
                        continue

                    # Переименовываем контиг
                    contig.id = f"{contig_basename}_{contig_name_parts[1]}"

                    # Формируем описание
                    if no_description:
                        contig.description = ""
                    else:
                        length = len(contig.seq)
                        spades_cov = contig_name_parts[5]
                        contig.description = f"length={length};spades_cov={spades_cov};{add_to_description}".rstrip(';')

                    # Записываем отфильтрованный контиг
                    SeqIO.write(contig, o_file, "fasta")

            return 0

        except FileNotFoundError as e:
            logger_asmb.error(f"Файл не найден: {e}", exc_info=True)
            return 1
        except PermissionError as e:
            logger_asmb.error(f"Ошибка доступа к файлу: {e}", exc_info=True)
            return 1
        except Exception as e:
            logger_asmb.error(f"Неожиданная ошибка при обработке контигов: {e}", exc_info=True)
            return 1
    else:
        print("В рабочей папке не найдены FASTA-файлы")
        exit_code = 127
    return exit_code


if __name__ == "__main__":
    min_len = sys.argv[1]
    contig_basename = sys.argv[2]
    no_description = sys.argv[3]
    add_to_description = sys.argv[4]
    exit(process_contigs())
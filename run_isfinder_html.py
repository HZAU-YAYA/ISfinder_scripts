#!/usr/bin/env python
# -*- coding: utf-8 -*-

import requests
from pyquery import PyQuery as pq
from bs4 import BeautifulSoup
import os
import logging
import argparse
import sys


LOG = logging.getLogger(__name__)
__version__ = "1.0.0"    #设置版本信息
__author__ = ("Boya Xu",)   #输入作者信息
__email__ = "834786312@qq.com"
__all__ = []


def add_help_args(parser):   #帮助函数
    parser.add_argument('--html', type=str, default=False, help="输入文件")
    parser.add_argument('--path', type=str, default=False, help="路径")
    parser.add_argument('--prefix', type=str, default=False, help="name")
    return parser

requests.packages.urllib3.disable_warnings()

def run_IS605_txt(html, prefix, path=os.getcwd()):
    f = open(html, 'r', encoding='utf-8')
    soup = BeautifulSoup(f, 'html.parser')
    result_table = soup.find('table', class_='result')
    header_row = result_table.find('tr')
    header_info = []
    header_cells = header_row.find_all('th')
    for cell in header_cells:
        header_info.append(cell.get_text())
    data_rows = result_table.find_all('tr')[1:]
    table_data = []
    for row in data_rows:
        row_cells = row.find_all('td')
        if row_cells:
            row_data = [cell.get_text() for cell in row_cells]
            links = row.find_all('a', target='_blank')
            link1 = links[0]['href'] if len(links) > 0 else ''
            link2 = links[1]['href'] if len(links) > 1 else ''
            link3 = links[2]['href'] if len(links) > 2 else ''
            # 添加链接到每组数据
            row_data.extend([link1, link2, link3])
            table_data.append(row_data)
    print(table_data)
    with open(path+'/'+prefix+'_out.tsv', 'w', encoding='utf-8') as tsvfile:
        # 写入表头
        tsvfile.write('\t'.join(header_info) + '\n')
        # 写入表格数据
        for row in table_data:
            tsvfile.write('\t'.join(row) + '\n')

def run_link(prefix, path=os.getcwd()):
    requests.packages.urllib3.disable_warnings()
    data = {}
    with open(path+'/'+prefix+'_out.tsv', 'r', encoding='utf-8') as file:
        next(file)
        lines = file.readlines()
    extracted_data = []
    for line in lines:
        columns = line.strip().split('\t')
        name = columns[1]
        url1 = columns[12]
        url2 = columns[13]
        url3 = columns[14]
        extracted_data.append([name, url1, url2, url3])
    print(extracted_data)
    for i in extracted_data:
        IS_name = i[0]
        url = i[1]
        try:
            response = requests.get(url, verify=False)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data: {str(e)}")
        soup = BeautifulSoup(response.text, 'html.parser')
        article = soup.find('article')
        host = article.find('div',class_='ascenseurAuto').text.strip()
        LE = article.find_all('span', class_='seq')[0].text.strip() if len(article) >0 else ''
        RE = article.find_all('span', class_='seq')[1].text.strip() if len(article) >0 else ''
        DNA = article.find_all('div', class_='seq')[0].text.strip() if len(article) >0 else ''
        for section in soup.find_all('section'):
            # 查找包含 ORF 编号的 <span> 元素
            orf = section.find('span', class_='entete_propriete', string='ORF number : ')
            if orf:
                # 获取 ORF 编号的值
                orf_n = orf.find_next_sibling(string=True).strip()
                if orf_n:
                    orf = section
                    break
        # orf = article.find_all('section')[5]
        # orf_n = orf.find_all('span', class_='entete_propriete', string='ORF number : ')[0].find_next_sibling(string=True).strip()
        data[IS_name] = [host, LE, RE, DNA, orf_n]
        if int(orf_n) > 0:
            for i in range(1, int(orf_n) + 1):
                oo = orf.find_all('span', class_='entete_propriete', string='ORF ')[i - 1]
                num = oo.find_next_sibling(string=True).strip()
                ORF_function = oo.find_next('span', class_='entete_propriete',string='ORF function : ').find_next_sibling(string=True).strip()
                na = oo.find_next('span', class_='entete_propriete', string='ORF function : ').find_next('span', class_='entete_propriete').find_next_sibling(string=True).strip()
                table = oo.find_next('table')
                rows = table.find_all('tr')
                th = rows[0].find_all('th')
                if len(th) != 5:
                    print('Error: The table format is non-standard')
                if th[0].text.strip() == 'Length':
                     length = rows[1].find_all('td')[0].text.strip() + '|' + rows[1].find_all('td')[1].text.strip()
                else:
                    length = ''
                    print('Error: length is not found')
                if th[1].text.strip() == 'Begin':
                     begin = rows[1].find_all('td')[2].text.strip()
                else:
                    begin = ''
                    print('Error: begin is not found')
                if th[2].text.strip() == 'End':
                     end = rows[1].find_all('td')[3].text.strip()
                else:
                    end = ''
                    print('Error: end is not found')
                if th[3].text.strip() == 'Strand':
                     strand = rows[1].find_all('td')[4].text.strip()
                else:
                    strand = ''
                    print('Error: strand is not found')
                if th[4].text.strip() == 'Fusion ORF':
                     fusion_orf = rows[1].find_all('td')[5].text.strip()
                else:
                    fusion_orf = ''
                    print('Error: fusion_orf is not found')
                if oo.find_next('p', class_='entete_propriete',string='ORF sequence : '):
                    ORF_protein = oo.find_next('div', class_='seq').text.strip()
                data[IS_name].append([num, ORF_function, na, ORF_protein, length, begin, end, strand, fusion_orf])
        elif int(orf_n) == 0:
            print(IS_name+' has no ORF')
        #pwd = os.getcwd()
        for II in data:
            folder_path = str(path+ '/' + prefix + '_db/' + II)
            os.makedirs(folder_path, exist_ok=True)
            if os.path.exists(folder_path):
                if int(data[II][4]) == 0:
                    with open(folder_path + '/' + II + '_seq.fasta', 'w') as fasta:
                        fasta.write('>' + II + '\n')
                        fasta.write(data[II][3])
                    with open(folder_path + '/' + II + '_LE_RE.tsv', 'w') as tsv:
                        tsv.write(II + '\n')
                        tsv.write('\n'.join(data[II][0:3]))
                elif int(data[II][4]) > 0:
                    with open(folder_path + '/' + II + '_seq.fasta', 'w') as fasta:
                        fasta.write('>' + II + '\n')
                        fasta.write(data[II][3])
                    with open(folder_path + '/' + II + '_LE_RE.tsv', 'w') as tsv:
                        tsv.write(II + '\n')
                        tsv.write('\n'.join(data[II][0:3]))
                    with open(folder_path + '/' + II + '_ORF' + '.tsv', 'w') as tsv:
                        tsv.write('IS\tGene\tORF_function\tlength\tStart\tEnd\tLength\tCoverage\n')
                    for iii in range(int(data[II][4])):
                        with open(folder_path + '/' + II + '_' + data[II][iii+5][2] + '_ORF' + str(iii + 1) + '.fasta', 'w') as fa:
                            fa.write('>' + II + '\t' + data[II][iii+5][2] + '\t' + data[II][iii+5][1] + '\n')
                            fa.write(data[II][iii+5][3])
                        with open(folder_path + '/' + II + '_ORF' +'.tsv', 'a') as tsv:
                            tsv.write(II + '\t' + data[II][iii+5][2] + '\t' + data[II][iii+5][1] + '\t' +
                                      data[II][iii+5][4] + '\t' + data[II][iii+5][5] + '\t' +
                                      data[II][iii+5][6] + '\t' + data[II][iii+5][7] + '\t' + data[II][iii+5][8]+'\n')
            else:
                print('ERROR:无法创建目录')
                break


def main():   #主函数，执行函数
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=''' 
name:statistic.py 
attention: python read_ISfinder --file  --pre  -p -l
version: %s
contact: %s <%s>\ 
静态网站爬取脚本，用于获取ISfinder的html脚本
''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()
    run_IS605_txt(args.html, args.prefix, args.path)
    run_link(args.prefix, args.path)

if __name__ == "__main__":           #固定格式，使 import 到其他的 python 脚本中被调用（模块重用）执行
    main()


# 搜索excel里的GSE，并搜索GEO页面
import re
import os
from openpyxl import Workbook
from pathlib import Path
import numpy as np
import pandas as pd
import urllib.response
from bs4 import BeautifulSoup
import requests
from tqdm import tqdm
import requests
from tqdm import tqdm
import time


# ----------------------------------------------------------------------------
# 从文件中搜索包含特定字符的行数据
def find_lines_in_file(filename, search_chars):
    matching_lines = []
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            for line_num, line in enumerate(file, 1):
                if search_chars in line:
                    matching_lines.append((line_num, line.strip()))
    except FileNotFoundError:
        print(f"文件 {filename} 未找到")
    return matching_lines

# 在excel中搜索包含特定字符的行数据
def find_rows_with_text(file_path, search_text, sheet_name=0, case_sensitive=False):
    """
    在 Excel 文件中查找包含特定文本的行

    参数:
    file_path: Excel 文件路径
    search_text: 要搜索的文本
    sheet_name: 工作表名称或索引，默认为第一个工作表
    case_sensitive: 是否区分大小写，默认为 False
    """
    # 读取 Excel 文件
    df = pd.read_excel(file_path, sheet_name=sheet_name)

    # 创建掩码，标识包含搜索文本的行
    if case_sensitive:
        mask = df.apply(lambda row: row.astype(str).str.contains(search_text).any(), axis=1)
    else:
        mask = df.apply(lambda row: row.astype(str).str.lower().str.contains(search_text.lower()).any(), axis=1)

    # 返回匹配的行
    result = df[mask]
    return result

# 在excel的特定列中搜索包含特定字符的行数据
def find_rows_in_column(file_path, search_text, column_name, sheet_name=0):
    """
    在 Excel 文件的特定列中查找包含特定文本的行
    """
    # 读取 Excel 文件
    df = pd.read_excel(file_path, sheet_name=sheet_name)

    # 检查列是否存在
    if column_name not in df.columns:
        print(f"列 '{column_name}' 不存在于工作表中")
        return pd.DataFrame()

    # 查找包含搜索文本的行
    mask = df[column_name].astype(str).str.contains(search_text, case=False)
    result = df[mask]

    return result


def append_to_excel_pandas(data, filename="data.xlsx", sheet_name="Sheet1"):
    """
    使用 pandas 追加数据到 Excel 文件

    参数:
    - data: 要追加的数据，可以是字典、列表或 DataFrame
    - filename: Excel 文件名
    - sheet_name: 工作表名
    """
    # 如果文件不存在，创建新文件
    if not os.path.exists(filename):
        df = pd.DataFrame(data)
        df.to_excel(filename, sheet_name=sheet_name, index=False)
        print(f"创建新文件 {filename} 并写入数据")
    else:
        # 读取现有文件
        existing_df = pd.read_excel(filename, sheet_name=sheet_name)

        # 将新数据转换为 DataFrame
        if isinstance(data, dict):
            new_df = pd.DataFrame([data])
        elif isinstance(data, list):
            new_df = pd.DataFrame(data)
        else:
            new_df = data

        # 合并数据
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)

        # 写回文件
        combined_df.to_excel(filename, sheet_name=sheet_name, index=False)
        # print(f"数据已追加到 {filename}")


def find_GSE_from_GEO(GSE):
    time.sleep(1)  # 延迟1秒

    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + str(GSE)
    response = requests.get(url)
    html_content = response.text
    data_GSE = GSE


    # Status
    pattern_Status = re.compile('<td>Status.*?<td>(.*?)</td>', re.S)
    try:
        data_Status = re.findall(pattern_Status, html_content)
    except IOError:
        data_Status = ['']
        print(f'ERROR：没有查询到 {GSE} 的Status')

    # Title
    pattern_Title = re.compile('<td nowrap>Title</td>.*?>(.*?)</td>', re.S)
    try:
        data_Title = re.findall(pattern_Title, html_content)
    except IOError:
        data_Title = ['']
        print(f'ERROR：没有查询到 {GSE} 的Title')

    # Organism
    pattern_Organism = re.compile('<td nowrap>Organism</td>\n<td>.*?>(.*?)</a>', re.S)
    try:
        data_Organism = re.findall(pattern_Organism, html_content)
    except IOError:
        data_Organism = ['']
        print(f'ERROR：没有查询到 {GSE} 的Organism')

    # Experiment type
    pattern_Experiment = re.compile('<td nowrap>Experiment type</td>\n<td>(.*?)<br>', re.S)
    try:
        data_Experiment = re.findall(pattern_Experiment, html_content)
    except IOError:
        data_Experiment = ['']
        print(f'ERROR：没有查询到 {GSE} 的Experiment')

    # summary
    pattern_Summary = re.compile('<td nowrap>Summary</td>\n<td style="text-align: justify">(.*?)<br>', re.S)
    try:
        data_Summary = re.findall(pattern_Summary, html_content)
    except IOError:
        data_Summary = ['']
        print(f'ERROR：没有查询到 {GSE} 的Summary')


    # design
    pattern_design = re.compile('<td nowrap>Overall design</td>\n<td style="text-align: justify">(.*?)<br>', re.S)
    try:
        data_design = re.findall(pattern_design, html_content)
    except IOError:
        data_design = ['']
        print(f'ERROR：没有查询到 {GSE} 的design')


    # Citation
    pattern_Citation = re.compile('<span class="pubmed_id" id="(.*?)">', re.S)
    try:
        data_Citation = re.findall(pattern_Citation, html_content)
        for num, item in enumerate(data_Citation):
            data_Citation[num] = '['+data_Citation[num]+']'+ r'(https://pubmed.ncbi.nlm.nih.gov/' + data_Citation[num] + r')'

    except IOError:
        data_Citation = ['NULL']
        print(f'ERROR：没有查询到 {GSE} 的Citation')


    # Email
    pattern_Email = re.compile('E-mail\(s\)</td>\n<td><a href="mailto:(.*?)">', re.S)
    try:
        data_Email = re.findall(pattern_Email, html_content)
    except IOError:
        data_Email = ['NULL']
        print(f'ERROR：没有查询到 {GSE} 的Email')


    # Platforms
    pattern_Platforms = re.compile('<a href="/geo/query/acc\.cgi\?acc=GPL(.*?)"', re.S)
    try:
        data_Platforms = re.findall(pattern_Platforms, html_content)
        for num, item in enumerate(data_Platforms):
            data_Platforms[num] = 'GPL'+data_Platforms[num]

    except IOError:
        data_Platforms = ['NULL']
        print(f'ERROR：没有查询到 {GSE} 的Platforms')


    # Samples
    pattern_Samples = re.compile('<a href="/geo/query/acc\.cgi\?acc=GSM(.*?)"', re.S)
    try:
        data_Samples = re.findall(pattern_Samples, html_content)
        for num, item in enumerate(data_Samples):
            data_Samples[num] = 'GSM' + data_Samples[num]
        # print(data_Samples_temp)
    except IOError:
        data_Samples = ['NULL']
        print(f'ERROR：没有查询到 {GSE} Samples')


    # BioProject
    pattern_BioProject = re.compile('<a href="https://www\.ncbi\.nlm\.nih\.gov/bioproject/.*?">(.*?)</a>', re.S)
    try:
        data_BioProject = re.findall(pattern_BioProject, html_content)
    except IOError:
        data_BioProject = ['NULL']
        print(f'ERROR：没有查询到 {GSE} 的BioProject')


    # Download family
    pattern_Download = re.compile(' target="_blank">(.*?)</a>', re.S)
    url_Download = re.compile('<a href="ftp://ftp\.ncbi\.nlm\.nih\.gov/geo/series/(.*?)" target="_blank">', re.S)
    try:
        data_Download = re.findall(pattern_Download, html_content)
        data_Download_url = re.findall(url_Download, html_content)
        # print(data_Download)
        # print(data_Download_url)
        for num, item in enumerate(data_Download):
            data_Download_url[num] = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/' + data_Download_url[num]
            data_Download[num] = '['+data_Download[num]+']'+ r'(' + data_Download_url[num] + r')'
    except IndexError:
        data_Download = ['NULL']
        print(f'ERROR：没有查询到 {GSE} 的Download')


    # Supplementary
    soup = BeautifulSoup(html_content, 'html.parser')
    temp = soup.find_all('table')  # 查找所有<h1>标签
    data_temp = temp[-2]
    # print(data_temp)
    # print(type(data_temp))
    try:
        pattern_Supplementary = re.compile('>GSE(.*?)</td>', re.S)
        data_Supplementary = re.findall(pattern_Supplementary, str(data_temp))
        pattern_Supplementary_url = re.compile('">'+GSE+'.*?</td>\n<td.*?">(.*?)</td>', re.S)
        data_Supplementary_url = re.findall(pattern_Supplementary_url, str(data_temp))
        # print(data_Supplementary[0])
        # print(data_Supplementary_url)
        for num, item in enumerate(data_Supplementary):
            data_Supplementary[num] = '[' + data_Supplementary[num] + ']' + r'(' + data_Supplementary_url[num] + r')'
    except IOError:
        data_Supplementary = ['NULL']
        print(f'ERROR：没有查询到 {GSE} 的design')


    all_data = [[data_GSE, data_Status, data_Title, data_Organism, data_Experiment, data_Summary, data_design, data_Citation, data_Email, data_Platforms, data_Samples, data_BioProject, data_Download, data_Supplementary]]

    return all_data
    # output_file = output_file_path + output_file_name
    #
    # append_to_excel_pandas(all_data,output_file, sheet_name="Sheet1")




# ----------------------------------------------------------------------------

print('GEO_research.py is running, let\'s go !!!')

input_file_path = './GEO_database/'
input_file_name = 'GEO_ALS.xlsx'
input_file = input_file_path + input_file_name
output_file_path = input_file_path
output_file_name = 'detail_'+input_file_name

try:
    df = pd.read_excel(input_file, sheet_name='Sheet1')
    print(r'Complete: {} file read'.format(input_file_name))
except IOError:
    print(r'ERROR: {} file not read'.format(input_file_name))
    df = 'DATA ERROR'


# 获取GSE数据
search_text = 'GSE'
search_df = find_rows_with_text(input_file, search_text)
if not search_df.empty:
    print(f"\nComplete：在 {input_file_name} 中找到 {len(search_df)} 行包含 {search_text} 的数据")
else:
    print(f"\nERROR: 没有在 {input_file_name} 中找到包含 '{search_text}' 的数据")

data_list = []
for line_num, line in search_df.iterrows():
    cell_data = search_df.at[line_num, 'GSE']  #指定列名的单元格
    data_list_temp = find_GSE_from_GEO(cell_data)
    data_list.extend(data_list_temp)
    print(f'Complete： {cell_data} 读写已完成，进度为：{line_num}/{len(search_df)-1}')

append_to_excel_pandas(data_list,output_file_path+output_file_name, sheet_name="Sheet1")
print(f'\n\nComplete： {output_file_name} 文件写入已完成')
# cell_data = search_df.at[0, 'GSE']
# find_GSE_from_GEO('GSE280753',output_file_path,output_file_name)





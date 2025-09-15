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

def print_mission(name):
    # 在下面的代码行中使用断点来调试脚本。
    print(f'Compete:  {name}')  # 按 Ctrl+F8 切换断点。

def test_file_open(file_path):
    """
    尝试以独占模式打开文件来判断是否被其他进程打开
    返回True表示文件可能被其他进程占用，False表示文件可用
    """
    try:
        # 尝试以追加模式打开文件，如果文件被其他进程占用会抛出异常
        with open(file_path, 'a', encoding='utf-8'):
            print(f"txt可以打开")
            pass
        return False
    except IOError:
        print(f"txt无法打开")
        return True


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

def GEO_spider(num,PM_txt,file_name):

    # cite
    pattern_cite = re.compile('<span class="position-number">(\d*?)</span>', re.S)
    data_cite = re.findall(pattern_cite, PM_txt)
    # print(data_cite)
    # print(f"num of cite: {len(data_cite)}")

    # Pubmed_ID
    pattern_PMID = re.compile('<span class="citation-part" style="display: none;">PMID: <span class="docsum-pmid">(\d*?)</span>', re.S)
    data_PMID = re.findall(pattern_PMID, PM_txt)
    # print(data_PMID)
    # print(f"num of PMID: {len(data_PMID)}")


    # DOI
    pattern_DOI = re.compile('<div class="Scholarscope_DOI">DOI: (.*?)</div>',
                               re.S)
    data_DOI = re.findall(pattern_DOI, PM_txt)
    try:
        data_DOI = 'https://doi.org/' + data_DOI[0]
    except IndexError:
        data_DOI = 'Not Found'
        # print(data_DOI)
    # print(f"num of DOI: {len(data_DOI)}")

    #Journal
    pattern_Journal = re.compile(
        '<div class="Scholarscope_JournalFrame notranslate"><div class="Scholarscope_Journal notranslate".*?;">(.*?)</div>',
        re.S)
    data_Journal = re.findall(pattern_Journal, PM_txt)
    # print(data_Journal)
    # print(f"num of Journal: {len(data_Journal)}")

    # IF
    pattern_IF = re.compile(
        '<div class="Scholarscope_Factor notranslate" style="background-color: .*?;">(.*?)</div>',
        re.S)
    data_IF = re.findall(pattern_IF, PM_txt)
    # print(data_IF)
    # print(f"num of IF: {len(data_IF)}")

    # CASR 中科院分区
    pattern_CASR = re.compile(
        '<div class="Scholarscope_Quartile notranslate" style="background-color:.*?;">(.*?)</div>',
        re.S)
    data_CASR = re.findall(pattern_CASR, PM_txt)
    # print(data_CASR)
    # print(f"num of CASR: {len(data_CASR)}")

    # YEAR
    pattern_YEAR = re.compile(
        '<div class="Scholarscope_Year">(.*?)</div>',
        re.S)
    data_YEAR  = re.findall(pattern_YEAR , PM_txt)
    # print(data_YEAR )
    # print(f"num of YEAR : {len(data_YEAR)}")


    # ArticleType
    pattern_Type = re.compile('<div class="Scholarscope_ArticleType notranslate">(.*?)</div>', re.S)
    data_Type = ['']
    try:
        data_Type = re.findall(pattern_Type, PM_txt)
    except IndexError:
        data_Type = ['']
    # print(data_Abstract)
    # print(f"num of Abstract: {len(data_Abstract)}")


    url = 'https://pubmed.ncbi.nlm.nih.gov/' + str(data_PMID[0]) + '/'
    # print(url)
    try:
        response = requests.get(url)
        html_content = response.text
        # print(html_content)
        # Title
        pattern_Title = re.compile('<h1 class="heading-title">(.*?)</h1>', re.S)
        try:
            data_Title = re.findall(pattern_Title, html_content)[0]
            data_Title.replace(r'\n', '')
            data_Title = [data_Title]
        except IOError:
            data_Title = ['']

        # Abstract
        pattern_Abstract = re.compile('<div class="abstract-content selected".*?<p>(.*?)</div>', re.S)
        try:
            data_Abstract = re.findall(pattern_Abstract, html_content)
        except IOError:
            data_Abstract = ['']
        # print(data_Abstract)
        # print(f"num of Abstract: {len(data_Abstract)}")
    except requests.exceptions.ConnectTimeout:
        data_Title = ['']
        data_Abstract = ['']
    except IndexError:
        data_Title = ['']
        data_Abstract = ['']
    except requests.exceptions.ReadTimeout:
        data_Title = ['']
        data_Abstract = ['']
    except requests.exceptions.RequestException as e:
        data_Title = ['']
        data_Abstract = ['']


    data_NUM = num
    # save data
    data = [[data_NUM, data_Type, data_Title,data_Abstract,data_Journal,data_CASR, data_IF, data_YEAR,  data_PMID, data_DOI, data_cite]]
    append_to_excel_pandas(data, "D:\\Research_Paper\\Paper\\1AD_MLVs\\2AD_mLVs_bioinformatics\\Pubmed_spider\\Pubmed_DATA\\"+file_name+".xlsx", sheet_name="Sheet1")


# 按装订区域中的绿色按钮以运行脚本。
if __name__ == '__main__':
    #url = paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GEO_id, sep="")
    #url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={} sep=\"\" ".format(GEO_id)
    #print_mission(url)
    # p = re.compile('GSE(\d*?)"', re.S)
    # print(p)
    file_name = "PIEZO1"
    file_path = "D:\\Research_Paper\\Paper\\1AD_MLVs\\2AD_mLVs_bioinformatics\\Pubmed_spider\\Pubmed_DATA\\" + file_name + ".txt"

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            txt = file.read()
            print(f"txt读取完成")
    except IOError:
        print(f"txt无法打开")

    # chipset = re.findall(p,txt)
    # print(chipset)
    # chipset = list(set(chipset))    # 去重复
    # chipset.sort()  # 排序
    # print(f"num of chipset: {len(chipset)}")

    # for j in chipset:
    #     print('正在爬取：GSE' + j)
    #     GEO_spider('GSE' + j, txt)

    # LIST
    # pattern_LIST = re.compile('<sp(.*?)type="checkbox" id=.*?>', re.S)
    pattern_LIST = re.compile('<article class="full-docsum"(.*?)</article>', re.S)
    data_LIST = re.findall(pattern_LIST, txt)
    # print(data_LIST)
    print(f"num of LIST: {len(data_LIST)}")

    for j in tqdm(range(len(data_LIST)),'proceed'):
        # print('正在爬取：LIST' + str(j+1))
        GEO_spider(str(j+1), data_LIST[j],file_name)

    # GEO_spider('GSE1', txt)



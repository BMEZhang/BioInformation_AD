import re
import os
from openpyxl import Workbook
from pathlib import Path
import numpy as np
import pandas as pd
import urllib.response
from bs4 import BeautifulSoup

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

def GEO_spider(GEO_GSE,GEO_txt):
    fyText = []
    fyTitle = []
    # GEO_txt = str(GEO_txt)
    # print(GEO_txt)
    # serials
    pattern_NUM = re.compile('an>(.*?).</span>', re.S)
    # print(pattern_NUM)
    data_NUM = re.findall(pattern_NUM, GEO_txt)
    # print(data_NUM)
    # print(f"num of serials: {len(data_NUM)}")

    # GSE
    # pattern_Title = re.compile('<p class="title">.*?justify">(.*?)</a></p>', re.S)
    pattern_GSE = re.compile('Accession: </dt><dd>(.*?)</dd></dl>', re.S)
    # print(pattern_GSE)
    data_GSE = re.findall(pattern_GSE, GEO_txt)
    # print(data_GSE)
    # print(f"num of GSE: {len(data_GSE)}")

    # Title
    # pattern_Title = re.compile('<p class="title">.*?justify">(.*?)</a></p>', re.S)
    pattern_Title = re.compile('<p class="title"><a href=.*?">(.*?)</a></p>', re.S)
    # print(pattern_Title)
    data_Title = re.findall(pattern_Title, GEO_txt)
    # print(data_Title)
    # print(f"num of Title: {len(data_Title)}")
    # data = data[0]
    # fyText.append(data)
    # fyTitle.append('Title')

    # Summary
    # pattern_Summary = re.compile('Summary</td>\n<td.*?justify">(.*?)<br></td>', re.S)
    pattern_Summary = re.compile('\(Submitter supplied\) (.*?)<a', re.S)
    # print(pattern_Summary)
    data_Summary = re.findall(pattern_Summary, GEO_txt)
    # print(data_Summary)
    # print(f"num of Summary: {len(data_Summary)}")
    # data = data[0]
    # fyText.append(data)
    # fyTitle.append('Summary')

    # GPL
    pattern_GPL = re.compile('Platform.*?">(GPL.*?)</a>', re.S)
    # print(pattern_GPL)
    data_GPL = re.findall(pattern_GPL, GEO_txt)
    # print(data_GPL)
    # print(f"num of GPL: {len(data_GPL)}")

    # DOWNLOAD_DATA
    pattern_DATA = re.compile('Download data: (.*?)</a>', re.S)
    # print(pattern_DATA)
    data_DATA = re.findall(pattern_DATA, GEO_txt)
    # print(data_DATA)
    # print(f"num of DATA: {len(data_DATA)}")

    # SAMPLES_SIZE
    pattern_SAMPLES = re.compile('<dt><a href="/gds/.*?">(.*?) Samples</a>', re.S)
    # print(pattern_SAMPLES)
    data_SAMPLES = re.findall(pattern_SAMPLES, GEO_txt)
    # print(data_SAMPLES)
    # print(f"num of SAMPLES: {len(data_SAMPLES)}")

    # Type
    pattern_Type = re.compile('Type: </dt><dd class="lng_ln">(.*?)</dd></dl>', re.S)
    # print(pattern_Type)
    data_Type = re.findall(pattern_Type, GEO_txt)
    # print(data_Type)
    # print(f"num of Type: {len(data_Type)}")

    data = [[data_NUM,data_GSE,data_Title,data_Summary,data_GPL,data_DATA,data_SAMPLES,data_Type]]

    # 保存excel
    # work_path = Path("D:\\Research_Paper\\Paper\\1AD_MLVs\\2AD_mLVs_bioinformatics\\spider_py\\GEO_database\\")
    # wb = Workbook()
    # gse_file = GEO_GSE
    # wb.save(filename=work_path.joinpath('',gse_file))
    # wb.save(filename=work_path.joinpath('备份','备份'+gse_file))
    append_to_excel_pandas(data, "D:\\Research_Paper\\Paper\\1AD_MLVs\\2AD_mLVs_bioinformatics\\spider_py\\GEO_database\\GEO_ALS.xlsx", sheet_name="Sheet1")
    # file_excel =  "D:\\Research_Paper\\Paper\\1AD_MLVs\\2AD_mLVs_bioinformatics\\spider_py\\GEO_database\\GEO_AD.xlsx"
    # df = pd.DataFrame(data, columns=['NUM', 'GSE', 'Title', 'Summary','GPL','DATA','SAMPLES','Type'])
    # df.to_excel(file_excel, index=False)

# 按装订区域中的绿色按钮以运行脚本。
if __name__ == '__main__':
    GEO_id = 'GSE5281'
    #url = paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GEO_id, sep="")
    #url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={} sep=\"\" ".format(GEO_id)
    #print_mission(url)
    p = re.compile('GSE(\d*?)"', re.S)
    print(p)
    file_path = "D:\\Research_Paper\\Paper\\1AD_MLVs\\2AD_mLVs_bioinformatics\\spider_py\\GEO_database\\ALS2.txt"

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
    pattern_LIST = re.compile('<sp(.*?)type="checkbox" id=.*?>', re.S)
    data_LIST = re.findall(pattern_LIST, txt)
    print(f"num of LIST: {len(data_LIST)}")

    for j in range(len(data_LIST)):
        print('正在爬取：LIST' + str(j+1))
        GEO_spider('GSE' + str(j+1), data_LIST[j])

    # GEO_spider('GSE1', txt)



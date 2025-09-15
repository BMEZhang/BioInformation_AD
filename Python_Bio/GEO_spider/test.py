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


url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295890'
response = requests.get(url)
html_content = response.text
soup = BeautifulSoup(html_content, 'html.parser')
titles = soup.find_all('table')  # 查找所有<h1>标签
print(titles[-2])
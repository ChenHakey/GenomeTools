from selenium import webdriver
from bs4 import BeautifulSoup
import re
import time
import pandas as pd

if __name__=='__main__':
    browser=webdriver.PhantomJS(executable_path=path)
    url='https://david.ncifcrf.gov/tools.jsp'
    browser.get(url)
    time.sleep(3)
    browser.find_element_by_id('LISTBox').send_keys(genesList)
    browser.find_element_by_xpath("//option[@value='REFSEQ_MRNA']").click()
    browser.find_element_by_xpath("//input[@name='rbUploadType']").click()
    browser.find_element_by_xpath("//input[@value='Submit List']").click()
    browser.find_elements_by_xpath("//li")[-2].click()
    browser.find_element_by_xpath("//button[@value='Functional Annotation Chart']").click()
    browser.switch_to_window(browser.window_handles[-1])
    browser.save_screenshot(r'C:\Users\Zheng Peng\Desktop\test.png')
    bsObj=BeautifulSoup(browser.page_source,'lxml')
    table=bsObj.find('table',{'class':'dataTable'})
    head=table.find('thead')
    columns=[]
    for i in head.findAll('th'):
        columns.append(i.get_text())
    body=table.find('tbody')
    data=[]
    for i in body.findAll('tr'):
        temp=[]
        for j in i.findAll('td'):
            temp.append(j.get_text())
        data.append(temp)
    df=pd.DataFrame(data,columns=columns)

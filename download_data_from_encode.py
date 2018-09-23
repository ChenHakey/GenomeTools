import requests
import re
import sys
import getopt
import os

def get_informations(url):
    infos={}
    response=requests.get(url)
    print('\t'.join(['Status code:',str(response.status_code)]))
    for i in response.content.strip().split('\n'):
        temp=i.split('\t')
        file_name=temp[0]
        temp_dic={}
        for j in temp[1].split('; '):
            key,value=j.split('=')
            temp_dic[key]=value
        infos[file_name]=temp_dic
    return infos

if __name__ == '__main__':
    cell_line=''
    modification_type=''
    out_dir=''
    doc='''
    Usage: python download_data_from_encode.py -c cell_line -t modification_type
    '''
    try:  
        opts,args=getopt.getopt(sys.argv[1:],"hc:t:o:",["help","cell_line","modification_type"])  
    except getopt.GetoptError:
        sys.exit()
    for opt,arg in opts:
        if opt in ('-h','--help'):
            print(doc)
            sys.exit()
        elif opt in ('-c','--cell_line'):
            cell_line=arg
        elif opt in ('-t','--modification_type'):
            modification_type=arg
    out_dir=os.getcwd()+'/'
    print('\t'.join(['cell_line:',cell_line]))
    print('\t'.join(['modification_type:',modification_type]))
    print('\t'.join(['output_directory:',out_dir]))
    print('\t'.join(['treatment','EtOH_0.02pct']))
    base_url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/'
    url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/files.txt'
    infos=get_informations(url)
    urls_list={}
    for file_name in infos:
        if cell_line==infos[file_name]['cell'] and modification_type==infos[file_name]['antibody'] and infos[file_name]['treatment']=='EtOH_0.02pct':
            urls_list[file_name]={}
            urls_list[file_name]['url']=base_url+file_name
            urls_list[file_name]['cell']=infos[file_name]['cell']
            urls_list[file_name]['antibody']=infos[file_name]['antibody']
            urls_list[file_name]['md5sum']=infos[file_name]['md5sum']
            urls_list[file_name]['size']=infos[file_name]['size']
        elif cell_line==infos[file_name]['cell'] and infos[file_name]['antibody']=='Control' and infos[file_name]['treatment']=='EtOH_0.02pct':
            urls_list[file_name]={}
            urls_list[file_name]['url']=base_url+file_name
            urls_list[file_name]['cell']=infos[file_name]['cell']
            urls_list[file_name]['antibody']=infos[file_name]['antibody']
            urls_list[file_name]['md5sum']=infos[file_name]['md5sum']
            urls_list[file_name]['size']=infos[file_name]['size']
    print('\t'.join(['number of files to download:',str(len(urls_list.keys()))]))
    for file_name in urls_list:
        print('\t'.join(['downloading...',file_name]))
        temp_url=urls_list[file_name]['url']
        md5sum=os.popen('wget -c %s > %s.downloading.log 2>&1 && md5sum %s'%(temp_url,file_name,file_name)).read().split()[0]
        temp_md5sum=urls_list[file_name]['md5sum']
        if temp_md5sum==md5sum:
            print('\t'.join([file_name,temp_md5sum,md5sum,'check','OK']))
        else:
            print('\t'.join([file_name,temp_md5sum,md5sum,'check','Failed']))



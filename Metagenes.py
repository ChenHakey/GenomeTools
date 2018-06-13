import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib import gridspec
import copy

def bisection(arr,aim):
    start=0
    end=len(arr)
    while start<end:
        mid=(start+end)//2
        if arr[mid]>aim:
            end=mid
        else:
            start=mid+1
    return start

def findRegion(arr,left,right):
    return [bisection(arr,left-1),bisection(arr,right)]

def findInterval(arr,left,right):
    return bisection(arr,right)-bisection(arr,left-1)

def read_peaks_file(path,chromList):
    peaks={}
    for chrom in chromList:
        peaks[chrom]=[]
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            if temp[0] in peaks:
                peaks[temp[0]].append(map(int,temp[1:3]))
    return peaks

def get_summits_of_peaks(path,chromList):
    summits={}
    for chrom in chromList:
        summits[chrom]=[]
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            if temp[0] in summits:
                summits[temp[0]].append(int(temp[1])+int(temp[4]))
    return summits

def get_summits_of_enhancers(enhancers,summits_of_peaks,chromList):
    summits_of_enhancers={}
    for chrom in chromList:
        summits_of_enhancers[chrom]=[]
    for chrom in enhancers:
        for start,end in enhancers[chrom]:
            summits=summits_of_peaks[chrom]
            left,right=findRegion(summits,start,end)
            summits_of_enhancers[chrom].extend(summits[left:right])
    return summits_of_enhancers

def read_bed_file(path,chromList):
    reads={}
    num_of_reads=0
    for chrom in chromList:
        reads[chrom]=[]
    test=0
    with open(path) as infile:
        for line in infile:
            if (test+1)%1000000==0:
                print('reading',test,'reads')
            test+=1
            temp=line.strip().split('\t')
            if temp[0] in reads:
                reads[temp[0]].append(int(temp[1]))
                num_of_reads+=1
    return reads,num_of_reads

def read_refseq_extract_tss(path,chromList):
    TSSs={}
    for chrom in chromList:
        TSSs[chrom]=[]
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split()
            if temp[2] in TSSs:
                if temp[3]=='+':
                    TSSs[temp[2]].append([int(temp[4]),temp[3]])
                else:
                    TSSs[temp[2]].append([int(temp[5]),temp[3]])
    for i in TSSs:
        TSSs[i]=sorted(TSSs[i],key=lambda x:x[0])
    return TSSs


def create_bins(peaks,TSSs,chromList,length=2500,bin_length=200):
    bins={}
    for chrom in chromList:
        bins[chrom]=[]
    for chrom in peaks:
        TSSs_loc=[i[0] for i in TSSs[chrom]]
        TSSs_strand=[i[1] for i in TSSs[chrom]]
        for start,end in peaks[chrom]:
            l,r=findRegion(TSSs_loc,start,end)
            for index,tss in enumerate(TSSs_loc[l:r]):
                temp=[[i,j-1] for i,j in zip(range(tss-length,tss+length,bin_length),range(tss-length+bin_length,tss+length+bin_length,bin_length))]
                if TSSs_strand[index]=='+':
                    bins[chrom].append(temp)
                else:
                    bins[chrom].append(temp[::-1])
    return bins

def count_it(bins,reads,num_of_reads,bin_length=200):
    counts=[]
    for chrom in bins:
        for bins_list in bins[chrom]:
            temp=[]
            for start,end in bins_list:
                count=findInterval(reads[chrom],start,end)
                temp.append(count*1000000000/(float(bin_length)*num_of_reads))
            counts.append(temp)
    counts=np.array(counts)
    return counts

'''

temp:[[------>]]
            ------>
                    ------->
                              ------->
                                --->  
'''


def merge(regions):
    merged_regions={}
    for chrom in regions:
        sorted_regions=copy.deepcopy(sorted(regions[chrom]))
        temp=[]
        for region in sorted_regions:
            if temp and temp[-1][1]>=region[0]-1:
                temp[-1][1]=max(temp[-1][1],region[1])       
            else:
                temp.append(region)
        merged_regions[chrom]=temp
    return merged_regions


def create_exclude_region(TSSs,exclude=2500):
    exclude_regions={}
    for chrom in TSSs:
        for tss,strand in TSSs[chrom]:
            if chrom in exclude_regions:
                exclude_regions[chrom].append([tss-exclude,tss+exclude])
            else:
                exclude_regions[chrom]=[[tss-exclude,tss+exclude]]
    exclude_regions=merge(exclude_regions)
    return exclude_regions

'''
overlap
(1)
region1  ------>
region2     ------->
(2)
region1     ------>
region2   ---------->
(3)
region1       ------>
region2   ------>
(4)
region1  ---------->
region2   ------>
non-overlap
(5)
region1  ----->
region2         ------>
(6)
region1          ----->
region2  ------>
'''

def overlap(region1,region2):
    status=None
    if min(region1[1],region2[1])-max(region1[0],region2[0])>0:
        if region2[1]>=region1[1]:
            if region2[0]>=region1[0]:
                status=1
            else:
                status=2
        else:
            if region1[0]>=region2[0]:
                status=3
            else:
                status=4
        return [max(region1[0],region2[0]),min(region1[1],region2[1])],status
    else:
        if region2[0]>region1[1]:
            status=5
        elif region1[0]>region2[0]:
            status=6
        return -1,status

def exclude_promoter(peaks,exclude_regions,chromList):
    excluded_peaks={}
    for chrom in chromList:
        excluded_peaks[chrom]=[]
    for chrom in peaks:
        temp_peaks=peaks[chrom]
        temp_regions=exclude_regions[chrom]
        index_of_peaks=0
        index_of_regions=0
        if len(temp_peaks)==0:
            excluded_peaks[chrom]=[]
            continue
        current_peak=temp_peaks[index_of_peaks]
        current_region=temp_regions[index_of_regions]
        while True:
            temp,status=overlap(current_region,current_peak)
            if temp!=-1:
                if status==1:
                    current_peak=[temp[1],current_peak[-1]]
                    index_of_regions+=1
                    if index_of_regions>len(temp_regions)-1:
                        break
                    current_region=temp_regions[index_of_regions]
                elif status==2:
                    excluded_peaks[chrom].append([current_region[0],temp[0]])
                    current_peak=[temp[1],current_peak[-1]]
                    index_of_regions+=1
                    if index_of_regions>len(temp_regions)-1:
                        break
                    current_region=temp_regions[index_of_regions]
                elif status==3:
                    excluded_peaks[chrom].append([current_peak[0],temp[0]])
                    current_region=[temp[1],current_region[-1]]
                    index_of_peaks+=1
                    if index_of_peaks>len(temp_peaks)-1:
                        break
                    current_peak=temp_peaks[index_of_peaks]
                elif status==4:
                    current_region=[temp[1],current_region[-1]]
                    index_of_peaks+=1
                    if index_of_peaks>len(temp_peaks)-1:
                        break
                    current_peak=temp_peaks[index_of_peaks]
            else:
                if status==5:
                    index_of_regions+=1
                    if index_of_regions>len(temp_regions)-1:
                        break
                    current_region=temp_regions[index_of_regions]
                else:
                    excluded_peaks[chrom].append(current_peak)
                    index_of_peaks+=1
                    if index_of_peaks>len(temp_peaks)-1:
                        break
                    current_peak=temp_peaks[index_of_peaks]
    return excluded_peaks

def create_bin_of_enhancers(enhancers,chromList,length=1000,bin_length=100):
    bins={}
    for chrom in chromList:
        bins[chrom]=[]
    for chrom in enhancers:
        for enhancer in enhancers[chrom]:
            if enhancer[1]-enhancer[0]<1000:
                continue
            mid_point=int((enhancer[0]+enhancer[1])/2)
            temp=[[start,end-1] for start,end in zip(range(mid_point-length,mid_point+length,bin_length),range(mid_point-length+bin_length,mid_point+length+bin_length,bin_length))]
            bins[chrom].append(temp)
    return bins

def create_bin_of_summits(summits,chromList,length=1000,bin_length=100):
    bins={}
    for chrom in chromList:
        bins[chrom]=[]
    for chrom in summits:
        for summit in summits[chrom]:
            temp=[[start,end-1] for start,end in zip(range(summit-length,summit+length,bin_length),range(summit-length+bin_length,summit+length+bin_length,bin_length))]
            bins[chrom].append(temp)
    return bins

if __name__ == '__main__':
    params = { 'axes.labelsize': '32', 
              'xtick.labelsize': '32', 
              'ytick.labelsize': '26', 
              'lines.linewidth': '4', 
              'legend.fontsize': '40', 
               'figure.figsize': '26, 24'}
    pylab.rcParams.update(params)
    chromList=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
               'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
               'chr20','chr21','chr22','chrX','chrY']
    plt.switch_backend('agg')
    
    print('reading refseq files...extact tss...')
    refseq_path=''
    TSSs=read_refseq_extract_tss(refseq_path,chromList)
    peaks_path_of_tumor=''
    peaks_path_of_normal=''

    print('reading peaks file...')
    peaks_of_tumor=read_peaks_file(peaks_path_of_tumor,chromList)
    peaks_of_normal=read_peaks_file(peaks_path_of_normal,chromList)

    summits_of_peaks_tumor=get_summits_of_peaks(peaks_path_of_tumor,chromList)
    summits_of_peaks_normal=get_summits_of_peaks(peaks_path_of_normal,chromList)

    print('creating bins...')
    bins_of_tumor_tss=create_bins(peaks_of_tumor,TSSs,chromList,2500,50)
    bins_of_normal_tss=create_bins(peaks_of_normal,TSSs,chromList,2500,50)
    bed_path_of_tumor=''
    bed_path_of_normal=''

    print('reading bed file(sorted)...')
    reads_of_tumor,num_of_reads_tumor=read_bed_file(bed_path_of_tumor,chromList)
    reads_of_normal,num_of_reads_normal=read_bed_file(bed_path_of_normal,chromList)

    print('creating excluding regions(TSS+/-2.5kb)...')
    exclude_regions=create_exclude_region(TSSs)

    print('excluding promoters(creating enhancers)...')
    enhancers_of_tumor=exclude_promoter(peaks_of_tumor,exclude_regions,chromList)
    enhancers_of_normal=exclude_promoter(peaks_of_normal,exclude_regions,chromList)

    summits_of_enhancers_tumor=get_summits_of_enhancers(enhancers_of_tumor,summits_of_peaks_tumor,chromList)
    summits_of_enhancers_normal=get_summits_of_enhancers(enhancers_of_normal,summits_of_peaks_normal,chromList)

    bins_of_tumor_enhancers=create_bin_of_summits(summits_of_enhancers_tumor,chromList,500,20)
    bins_of_normal_enhancers=create_bin_of_summits(summits_of_enhancers_normal,chromList,500,20)

    print('counting reads...')
    counts_of_tumor_tss=count_it(bins_of_tumor_tss,reads_of_tumor,num_of_reads_tumor,50)
    counts_of_normal_tss=count_it(bins_of_normal_tss,reads_of_normal,num_of_reads_normal,50)
    counts_of_tumor_enhancers=count_it(bins_of_tumor_enhancers,reads_of_tumor,num_of_reads_tumor,20)
    counts_of_normal_enhancers=count_it(bins_of_normal_enhancers,reads_of_normal,num_of_reads_normal,20)

    print('drawing...')
    plt.figure(figsize=(40,40))
    gs=gridspec.GridSpec(4,4,height_ratios=[2,5,10,1],width_ratios=[1,1,1,1])
    
     
    TSSs_of_tumor=counts_of_tumor_tss.mean(axis=0)
    TSSs_of_normal=counts_of_normal_tss.mean(axis=0)
    ESs_of_tumor=counts_of_tumor_enhancers.mean(axis=0)
    ESs_of_normal=counts_of_normal_enhancers.mean(axis=0)


    index_of_tumor_tss=np.argsort(counts_of_tumor_tss.sum(axis=1))
    index_of_normal_tss=np.argsort(counts_of_normal_tss.sum(axis=1))
    index_of_tumor_enhancer=np.argsort(counts_of_tumor_enhancers.sum(axis=1))
    index_of_normal_enhancer=np.argsort(counts_of_normal_enhancers.sum(axis=1))


    plt.subplot(gs[0,:])
    plt.barh(0,10,1,0,color='black')
    title='H3K27ac profiles in Lung cancer cells & normal cells'
    plt.text(5,5,title,fontsize=50,ha='center',va='top',wrap=True)
    plt.ylim(0,10)
    plt.xlim(0,10)
    plt.axis('off')


    plt.subplot(gs[1,0])
    plt.title('Promoter of Tumor(ID:2573)',size=40)
    plt.plot([i for i in range(len(TSSs_of_tumor))],TSSs_of_tumor)
    plt.xticks([0,50,100],['-2.5K','TSS','+2.5k'],rotation=0)

    plt.subplot(gs[1,1])
    plt.title('Promoter of Normal(ID:2573)',size=40)
    plt.plot([i for i in range(len(TSSs_of_normal))],TSSs_of_normal)
    plt.xticks([0,50,100],['-2.5K','TSS','+2.5k'],rotation=0)

    plt.subplot(gs[1,2])
    plt.title('Enhancer of Tumor(ID:2573)',size=40)
    plt.plot([i for i in range(len(ESs_of_tumor))],ESs_of_tumor)
    plt.xticks([0,25,50],['500bp','summit','500bp'],rotation=0)

    plt.subplot(gs[1,3])
    plt.title('Enhancer of Normal(ID:2573)',size=40)
    plt.plot([i for i in range(len(ESs_of_normal))],ESs_of_normal)
    plt.xticks([0,25,50],['500bp','summit','500bp'],rotation=0)

    plt.subplot(gs[2,0])
    data=np.log10((counts_of_tumor_tss+1))[index_of_tumor_tss[::-1],:]
    plt.imshow(data, cmap='Reds',interpolation='nearest',aspect='auto',vmin=0,vmax=1.0)

    plt.subplot(gs[2,1])
    data=np.log10((counts_of_normal_tss+1))[index_of_normal_tss[::-1],:]
    plt.imshow(data, cmap='Reds',interpolation='nearest',aspect='auto',vmin=0,vmax=1.0)

    plt.subplot(gs[2,2])
    data=np.log10((counts_of_tumor_enhancers+1))[index_of_tumor_enhancer[::-1],:]
    plt.imshow(data, cmap='Reds',interpolation='nearest',aspect='auto',vmin=0,vmax=1.0)

    plt.subplot(gs[2,3])
    data=np.log10((counts_of_normal_enhancers+1))[index_of_normal_enhancer[::-1],:]
    plt.imshow(data, cmap='Reds',interpolation='nearest',aspect='auto',vmin=0,vmax=1.0)


    ax1=plt.subplot(gs[3,3])
    cmap=matplotlib.cm.Reds
    norm=matplotlib.colors.Normalize(vmin=0,vmax=1.0)
    cb1=matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal',ticks=[0,0.5,1])
    cb1.set_label('$log_{10}(RPKM$)')
    cb1.ax.set_xticklabels(['0','0.5','1'])

    print(num_of_reads_tumor,num_of_reads_normal)
    print('done')
    plt.savefig('test.pdf')
    plt.close('all')













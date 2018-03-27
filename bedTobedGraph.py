import sys, getopt


'''
chromosome sizes
'''
chromSize={
    'chr1':249250621,
    'chr2':243199373,
    'chr3':198022430,
    'chr4':191154276,
    'chr5':180915260,
    'chr6':171115067,
    'chr7':159138663,
    'chrX':155270560,
    'chr8':146364022,
    'chr9':141213431,
    'chr10':135534747,
    'chr11':135006516,
    'chr12':133851895,
    'chr13':115169878,
    'chr14':107349540,
    'chr15':102531392,
    'chr16':90354753,
    'chr17':81195210,
    'chr18':78077248,
    'chr20':63025520,
    'chrY':59373566,
    'chr19':59128983,
    'chr22':51304566,
    'chr21':48129895      
}
'''
index:[0,1,2,3,4, 5, 6, 7]
vlaue:[1,3,5,7,9,11,13,15]
              |
              8(4)
    time complexity:O(log(n))
'''
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
'''
index:[0,1,2,3,4, 5, 6, 7]
vlaue:[1,3,5,7,9,11,13,15]
      [1,3,5,7]+[8]+[9,11,13,15]
      memory copy(move element):O(n)
'''
def insertion(interval,depth,value):
    insertPos=bisection(interval,value)
    if interval[insertPos-1]==value:
        insertPos-=1
    else:
        interval.insert(insertPos,value)
        depth.insert(insertPos,depth[insertPos-1])
    return interval,depth,insertPos
'''
           2:4  
[1,2,3,0,1]==>[1,2,3+1,0+1,1]
'''    
def spread(depth,start,end):
    for i in range(start,end):
        depth[i]+=1
    return depth


def transform(interval,depth,chrom,scale=1):
    preItem=[]
    depth=[i*scale for i in depth]
    for i in range(len(depth)):
        if depth[i]==0:
            continue
        elif preItem and interval[i]==preItem[2] and depth[i]==preItem[3]:
            preItem=[chrom,preItem[1],interval[i+1],depth[i]]
        elif not preItem:
            preItem=[chrom,interval[i],interval[i+1],depth[i]]
        else:
            yield preItem
            preItem=[]
            preItem=[chrom,interval[i],interval[i+1],depth[i]]
    yield preItem

if __name__ == '__main__':
    inputPath=''
    try:
        opts,args=getopt.getopt(sys.argv[1:],':i:')
    except getopt.GetoptError:
        print('test.py -i <inputfile> > <outputfile>')
        sys.exit(2)
    for opt,arg in opts:
        if opt in ("-i", "--ifile"):
            inputPath=arg
    with open(inputPath) as inFile:
        scale=1
        preChrom=''
        depth=[]
        interval=[]
        for i in inFile:
            temp=i.strip().split()
            chrom=temp[0]
            start=int(temp[1])
            end=int(temp[2])               
            if preChrom=='':
                preChrom=chrom
                interval=[0,chromSize[chrom]]
                depth=[0]
                interval,depth,insertPoss=insertion(interval,depth,start)
                interval,depth,insertPose=insertion(interval,depth,end)
                depth=spread(depth,insertPoss,insertPose)
            elif preChrom==chrom:    
                interval,depth,insertPoss=insertion(interval,depth,start)
                interval,depth,insertPose=insertion(interval,depth,end)
                depth=spread(depth,insertPoss,insertPose)
            else:
                for i in transform(interval,depth,preChrom,scale):
                    print('\t'.join(map(str,i)))
                preChrom=chrom
                interval=[0,chromSize[chrom]]
                depth=[0]
        for i in transform(interval,depth,preChrom,scale):
            print('\t'.join(map(str,i))) 

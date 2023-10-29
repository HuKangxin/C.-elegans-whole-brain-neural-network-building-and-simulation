#Created by Kangxin Hu at Beihang University, on Oct.2023.

from statistics import mean
import random
import math
def decode1(the_list):
    action = 0
    action_pos=[0]
    thres = -35
    for i in range(0,len(the_list)-1):
        if the_list[i] < thres:
            if the_list[i+1] >= thres:
                action += 1
        if the_list[i] > thres:
            if the_list[i+1] <= thres:
                action_pos.append(i)
    result=0.01*action+mean(action_pos)/40000
    return result

def limit(x):
    if x>=0.1:
        x=0.3*x
    return x

def decode2(the_list):
    action=0
    thres=-50
    for i in range(0,len(the_list)-1):
        if the_list[i+1]>=thres:
            action=1
            break
    return action

def act_record(voltage):
    record=[]
    record2=[]
    pos_0=[]
    pos_1 = []
    pos_t=[]
    for i in range(0,150):
        record.append(0)
    for i in range(0,len(voltage)):
        record2.append(0)
    for j in range(0,len(voltage)):
        if voltage[j]>=-60:
            record2[j]=1
    if len(voltage) > 150:
        for i in range(0,len(record2)-1):
            if record2[i] == 0:
                if record2[i+1] == 1:
                    pos=math.floor((i+1)*150/len(record2))
                    pos_0.append(pos)
            if record2[i] == 1:
                if record2[i+1] == 0:
                    pos = math.floor((i + 1) * 150 / len(record2))
                    pos_1.append(pos)
        for i in range(0,len(pos_1)):
            pos_t.append(pos_0[i])
            pos_t.append(pos_1[i])
        pos_t.append(150)
        pos_t.insert(0,0)
        print(pos_t)
        for i in range(1,len(pos_t)):
            if i%2==1:
                for j in range(pos_t[i-1],pos_t[i]):
                    record[j]=0
            if i % 2 == 0:
                for j in range(pos_t[i-1],pos_t[i]):
                    record[j]=1
    return record
'''
def act_record2(voltage):
    record=[]
    split=[]
    count=0
    for i in range(0,len(voltage)):
        if voltage[i]<=-65:
            record[i]=0
        else:
            record[i]=1
    for j in range(0,len(record)-1):
        if record[j]==1:
            split.append(0)
            if record[j+1]==1:
                split[count]+=1
            elif record[j+1]==0:
                count+=2
                split.append(j)
    for k in range(0,len(split)):
        if k%2==0:
            for x in range(split[k+1]-split[k],split[k+1]+1):
                record[x]=(x-split[k+1]+split[k])/split[k]
    return record
'''
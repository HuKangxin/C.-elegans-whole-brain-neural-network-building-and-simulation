# Created by Kangxin Hu at Beihang University, on Oct.2023.

import struct
import numpy as np


def bytesToFloat(by):
    ba = bytearray(by)
    return struct.unpack("f", ba)[0]


def floatToBytes(f):
    bs = struct.pack("f", f)
    return list(bs)


def Hex(x):
    temp = hex(x)
    if len(temp) == 3:
        return "0%s" % temp[-1]
    return temp[-2:]


def f2str(fs):   # 浮点数数组到字符串
    result = ''
    for f in fs:
        temp = floatToBytes(f)
        for i in temp:
            result += Hex(i)
    return result


def str2f(s):  # 字符串到浮点数数组
    fs = []
    result = []
    numbers = int(len(s)/8)
    for i in range(numbers):
        fs.append(s[i*8:i*8+8])
    for f in fs:
        hexs = bytes.fromhex(f)
        result.append(bytesToFloat(hexs))
    return np.array(result)

def zhuanstr(s):
    result=''
    for f in s:
        f=round(f,4)
        x='{:.4f}'.format(f)
        num_0=20-len(x)
        for i in range(0,num_0):
            x=x+'0'
        result+=x
    return result

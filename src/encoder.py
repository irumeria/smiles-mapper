import numpy as np
import json
import pandas as pd
import re
from collections import deque
import copy

# 广度优先搜索--分形版
# 分子结构相当于一个有环无向图
# 使用迭代次数来解决循环问题
# 记录历史来防止往回走
# 最重要的是：把search_queue里面每个节点都和到达它的时候经过的历史节点数组绑定在一起,做出 path
# 因为这里搜索主链，要找的应该是最长的，因此必须等search_queue搜索结束
# 然后对比每一个能到达目标节点的path，设最长那一条为主链
def search_main_chain(nodeMap, startIndex,endIndex):
    search_queue = []
    init_path = []
    init_path.append(startIndex)
    search_queue.append(init_path)
    target_path = [] # 用来存放能够到达目标node的path
    while search_queue:
        path = search_queue.pop(0)
        if not path[-1] in path[:-1]: # 不能走回头路
            # print(path)
            if path[-1] == endIndex: # 找到了！
                target_path.append(path) 
            else: # 没找到，将这个节点的邻居小伙伴作为结尾的新path加入队列中
                for newNode in nodeMap[path[-1]]["neighbors"]:
                    newPath = copy.deepcopy(path)
                    newPath.append(newNode["index"])
                    search_queue.append(newPath)
        pass
    if len(target_path) > 0:
        lenPath = []
        for tp in target_path:
            lenPath.append(len(tp))
        maxLenIndex = lenPath.index(max(lenPath))
        mainChain = target_path[maxLenIndex] # 主链是所有之中最长的那个
        return mainChain
    else: # 没有到达目的地的路径
        return []      

# 迭代找出所有支链/支链的支莲/...的函数
# 支链有两个东西标志尽头，一个是只有一个邻原子，另外一个是回到 baseChain 上
# baseChain 一开始是主链，迭代一遍之后变成主链+次级支链，以此类推
# startIndex： 从baseChain上某个节点开始搜索支链
def search_chain_iterative(nodeMap, baseChain, startIndex):
    loopers = [] # 记录分子中的环
    chainMap = {}
    search_queue = []
    for neighbor in nodeMap[startIndex]["neighbors"]:
        search_queue.append([neighbor["index"]])
    target_path = [] # 用来存放能够到达目标node的path
    while search_queue:
        path = search_queue.pop(0)
        if not path[-1] in path[:-1]: # 不能走回头路
            # 找到了！
            if len(nodeMap[path[-1]]["neighbors"]) == 1 and path[-1] not in baseChain:
                target_path.append(path) 
            elif path[-1] in baseChain : 
                if path[-1] != startIndex and len(path) > 1:
                    path.pop()
                    target_path.append(path) 
            else: # 没找到，将这个节点的邻居小伙伴作为结尾的新path加入队列中
                for newNode in nodeMap[path[-1]]["neighbors"]:
                    newPath = copy.deepcopy(path)
                    newPath.append(newNode["index"])
                    search_queue.append(newPath)
        pass
    if len(target_path) > 0:
        lenPath = []
        for tp in target_path:
            lenPath.append(len(tp))
        maxLenIndex = lenPath.index(max(lenPath))
        mainChain = target_path[maxLenIndex] # 主链是所有之中最长的那个
        
        # 来找环！
        # 判定环的标准：
        # 1. 次级支链上的原子连接到了更初级的支链/主链上
        # 2. 同一支链内，出现跨index的连接
        for x in mainChain:
            xIndex = mainChain.index(x) 
            for neighbor in nodeMap[x]["neighbors"]:
                if neighbor["index"] in baseChain and neighbor["index"] != startIndex:
                    loop = []
                    loop.append(neighbor["index"]) # 把baseChain上面那个原子吐出来
                    loop.append(x)
                    loopers.append(loop)
                for y in mainChain:
                    yIndex = mainChain.index(y)
                    if y == neighbor["index"] and yIndex != xIndex - 1 and yIndex != xIndex+1:
                        loop = []
                        
                        contained_flag = False
                        for item in loopers:
                            if item[0] == x and item[1] == y:
                                contained_flag = True
                        if not contained_flag:
                            loop.append(y) # 把baseChain上面那个原子吐出来
                            loop.append(x)
                            loopers.append(loop)
        
        for x in mainChain:
            baseChain.append(x)
        
                
        for x in mainChain:
            chainMap[str(x)],baseChain,newLoopers = search_chain_iterative(nodeMap, baseChain, x)
        if len(newLoopers) > 0:
            for x in newLoopers:
                loopers.append(x)

    
    return chainMap,baseChain,loopers
    
# 利用 chainMap 记录的内容，递归地将所有支链整合进smiles字符数组内
def build_brunch(nodeMap,chainMap):
    smiles=[]
    smileIndexs=[] 
    for key in chainMap:
        index = int(key)
        if len(nodeMap[index]["atom"]) == 1:
            smiles.append(nodeMap[index]["atom"])
        else:
            smiles.append("["+nodeMap[index]["atom"]+"]")
        smileIndexs.append(index)
        if len(chainMap[key]) > 0: # 它底下还有更次级的支链！
            brunchMap = chainMap[key]
            brunchSmiles,burnchIndexs = build_brunch(nodeMap=nodeMap,chainMap=brunchMap)
            
            smiles.append("(")
            smileIndexs.append(-1) # 标志 "(" 
            for i in range(0,len(brunchSmiles)):
                if burnchIndexs[i] == -1:
                    smiles.append("(")
                elif burnchIndexs[i] == -2:
                    smiles.append(")")
                else:
                    smiles.append(brunchSmiles[i])
                smileIndexs.append(burnchIndexs[i])
            smiles.append(")")
            smileIndexs.append(-2) # 标志 ")"
    
    return smiles,smileIndexs
    

def map_to_smiles(nodeMap):
    mainChain = search_main_chain(nodeMap, 0, len(nodeMap)-1)
    # 插入支链
    chainMap = {}
    baseChain = copy.deepcopy(mainChain)
    loopers = []
    print(mainChain)
    for x in mainChain:
        xIndex = mainChain.index(x) 
        for neighbor in nodeMap[x]["neighbors"]:
            for y in mainChain:
                yIndex = mainChain.index(y)
                if y == neighbor["index"] and yIndex != xIndex - 1 and yIndex != xIndex+1:
                    loop = []
                    contained_flag = False
                    for item in loopers:
                        if item[0] == x and item[1] == y:
                            contained_flag = True
                    if not contained_flag:
                        loop.append(y) # 把baseChain上面那个原子吐出来
                        loop.append(x)
                        loopers.append(loop)

    for index in range(len(mainChain)-1,-1,-1):
        x = mainChain[index]
        chainMap[str(x)],baseChain,newLoopers = search_chain_iterative(nodeMap, baseChain, x)
        if len(newLoopers) > 0:
            for newLoop in newLoopers:
                loopers.append(newLoop)
    print(chainMap,loopers)
    
    chainMap =dict(reversed(list(chainMap.items())))
    # smileIndex: 记录打造的smiles里面符号在nodeMap中index, 与 smiles char数组下标一一对应
    smiles,smileIndex = build_brunch(nodeMap=nodeMap,chainMap=chainMap)
    print(smiles)
    
    # 标记环
    counter = 0
    for loop in loopers:    
        counter += 1    
        if loop[0] > loop[1]:
            loop.reverse()
        for i in range(0,len(smiles)):
            if smileIndex[i] == loop[0]:
                smiles.insert(i+1, str(counter))
                smileIndex.insert(i+1, -3)
                break
        for i in range(0,len(smiles)):
            if smileIndex[i] == loop[1]:
                smiles.insert(i+1, str(counter))
                smileIndex.insert(i+1, -3)
                break
    

    
    # 标记单键/双键
    i = 0
    while i < len(smiles):
        if smileIndex[i] > 0:
            for neighbor in nodeMap[smileIndex[i]]["neighbors"]:
                if neighbor["bond"] == 2 and neighbor["index"] < smileIndex[i]:
                    smiles.insert(i, "=")
                    smileIndex.insert(i, -4)
                    i += 1
                elif neighbor["bond"] == 3 and neighbor["index"] < smileIndex[i]:
                    smiles.insert(i, "#")
                    smileIndex.insert(i, -5)
                    i += 1
        i += 1

    smiles = "".join(smiles)

    print(smiles)   
    
    return smiles 






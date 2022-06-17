import numpy as np
import json
import pandas as pd
import re

# TODO:  1.支持手性符号的解析 \ & @


def smiles_to_map(smiles):

    smiles = np.array(list(smiles))  # string 转 numpy char array

    # 这个函数不太优雅.
    def piBondConveter(nodeMap):
        for node in nodeMap:
            for neighbor in node['neighbors']:
                
                neighbor['atom'] = nodeMap[int(
                    neighbor['index'])]['atom']  # 动点手脚，顺便加个atom属性
                if neighbor['atom'] == node['atom'] == 'c':
                    neighbor['bond'] = 4  # 更改为大PI键
                
                if neighbor['bond'] != 1: # 由于我奇怪的history机制，双键只加到了后面的那个原子上面...
                    for neigh in nodeMap[neighbor['index']]['neighbors']:
                        if neigh['index'] == node['index']:
                            neigh['bond'] = neighbor['bond']
                            
        return nodeMap

    def buildNode(smChar, index, preIndex, i, smiles, historyBox):

        node = {}
        node['atom'] = smChar
        node['index'] = index

        # 判断原子的邻居
        if not ifAtomBehind(smiles, i) and i == 0:
            node['neighbors'] = []
        elif not ifAtomBehind(smiles, i):
            # 到达smiles尾部
            node['neighbors'] = [
                {"index": preIndex, "bond": historyBox['bond'][0]}]
        elif i == 0:
            node['neighbors'] = [
                {"index": index+1, "bond": historyBox['bond'][1]}]
        else:
            node['neighbors'] = [{"index": preIndex, "bond": historyBox['bond'][0]}, {
                "index": index+1, "bond": historyBox['bond'][1]}]

        return node

    def ifAtomBehind(smiles, i):
        for smChar in smiles[i+1:]:
            if re.match(r"[A-Za-z]", smChar):
                return True
        return False

    def analyse(smiles, i=0, index=0, nodeMap=[], preIndex=-1, historyBox=None):

        smChar = smiles[i]

        # 保存环，分支等几何结构的关键点，以空间换时间
        if historyBox is None:
            historyBox = {}
            # 存放分支起始点在nodeMap的序号,检索时自动拿取最后一个位置的，用完要删去
            historyBox['branch'] = []
            # 存放各个循环头的序号，不去删除其中的元素，用数组下标表示smiles里面的 1 2 等数字
            historyBox['loops'] = []
            # 存放各个长元素头的 i , 注意不是index了，这是为了记录原子的字符长度
            historyBox['longAtom'] = None
            # 存放某个原子前后两个键数,0:左边，1：右边
            # 长度2足够解析线性分子编码了
            historyBox['bond'] = [1, 1]

        # 处理键
        # 新读取到的键从[1]进入，从[0]出去
        historyBox['bond'][0] = historyBox['bond'][1]
        if re.match(r"[=]", smChar):
            historyBox['bond'][1] = 2  # 双键
        elif re.match(r"[#]", smChar):
            historyBox['bond'][1] = 3  # 三键
        else:
            historyBox['bond'][1] = 1  # 单键

        # 添加特殊结构导致的的邻节点
        # preIndex : 指着前一个原子节点在nodeMap的序号，是这里遇到smiles特殊符号时，操作位置的依据
        # 抵达环结尾：目前只支持最多九个环（因为懒）
        if re.match(r"[1-9]", smChar) and len(historyBox["loops"]) >= int(smChar):
            loopHeader = historyBox["loops"][int(
                smChar)-1]  # 环开头原子节点序号应为数字所在序号
            # for x in nodeMap[preIndex]['neighbors']:
            # if x['index'] == index:
            # nodeMap[preIndex]['neighbors'].remove(x) # 删除原右节点
            nodeMap[preIndex]['neighbors'].append(
                {"index": loopHeader, "bond": historyBox['bond'][0]})
            nodeMap[loopHeader]['neighbors'].append(
                {"index": preIndex, "bond": historyBox['bond'][0]})
            # if ifAtomBehind(smiles=smiles,i=i):
            # nodeMap[loopHeader]['neighbors'].append({"index":index,"bond":historyBox['bond'][0]})
            # preIndex = loopHeader
        elif re.match(r"[)]", smChar):
            # 分支头 前面一个原子的 nodeMap序号 应为 "(" 的前一个
            # print(nodeMap[preIndex], index)
            for x in nodeMap[preIndex]['neighbors']:
                if x['index'] == index:
                    nodeMap[preIndex]['neighbors'].remove(x)
            preIndex = historyBox["branch"].pop() - 1
            nodeMap[preIndex]['neighbors'].append(
                {"index": index, "bond": historyBox['bond'][0]})
        elif re.match(r"[\]]", smChar):
            atomName = ''.join(smiles[historyBox['longAtom']+1:i])
            node = buildNode(smChar=atomName, index=index,
                             preIndex=preIndex, i=i, smiles=smiles, historyBox=historyBox)
            nodeMap.append(node)  # 这里nodeMap的下标必须严格和index - 1对应上
            index += 1
            preIndex = index - 1
            historyBox['longAtom'] = None

        if re.match(r"[(]", smChar):  # 记录分支
            historyBox['branch'].append(index)
        elif re.match(r"[1-9]", smChar):  # 记录环
            historyBox['loops'].append(index - 1)
        elif re.match(r"[\[]", smChar):
            historyBox['longAtom'] = i  # 记录长原子

        # 筛选原子节点的规则
        # 1 是个字母
        # 2 不在 [] 内， 在的话，之后遇到 ] 后在上文代码中一起添加
        if re.match(r"[A-Za-z]", smChar) == None or historyBox['longAtom'] != None:
            if i == len(smiles) - 1:
                return nodeMap
            else:
                # 相当于跳过这个节点，index不增加
                return analyse(smiles, i=i+1, index=index, preIndex=preIndex, nodeMap=nodeMap, historyBox=historyBox)

        node = buildNode(smChar=smChar, index=index,
                         preIndex=preIndex, i=i, smiles=smiles, historyBox=historyBox)
        nodeMap.append(node)  # 这里nodeMap的下标必须严格和index - 1对应上

        # 必须遇到原子节点，这两个才会上升
        index += 1
        preIndex = index - 1

        if i == len(smiles) - 1:
            return nodeMap
        else:
            return analyse(smiles, i=i+1, index=index, preIndex=preIndex, nodeMap=nodeMap, historyBox=historyBox)

    nodeMap = analyse(smiles)
    return piBondConveter(nodeMap)

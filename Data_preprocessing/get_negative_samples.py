import random


def get_PP_noninteraction():
    with open("9606.protein.links.v11.0.txt", "r") as f:  # 打开文件
        PP_interaction = f.read()  # 读取PPI文件
        # print(PP_interaction)
    PP_interaction_lines = PP_interaction.split("\n")
    # PP_interaction_lines_list=PP_interaction_lines[0].split(" ")
    # print(PP_interaction_lines_list)

    # all_p=[]
    # for i in PP_interaction_lines:
    #     PP_interaction_lines_list=i.split(" ")
    #     if(len(PP_interaction_lines_list)>1):
    #         if(PP_interaction_lines_list[0] in all_p):
    #             pass
    #         else:
    #             all_p.append(PP_interaction_lines_list[0])
    #         if(PP_interaction_lines_list[1] in all_p):
    #             pass
    #         else:all_p.append(PP_interaction_lines_list[1])
    # print("all_p的长度为：")
    # print(len(all_p))

    with open("all-protein_id.txt", "r") as f:  # 打开文件
        protein_id = f.read()  # 读取筛选后所得的所有蛋白质文件
        # print(PP_interaction)
    all_p = protein_id.split("\n")

    PM = [[0] * len(all_p) for i in range(len(all_p))] # Filling adjacency matrix with 0

    print("二维数组的shape为")
    print(len(PM))
    print(len(PM[0]))

    PP_interaction=[]
    for i in PP_interaction_lines:
        PP_interaction_lines_list = i.split(" ")
        if (len(PP_interaction_lines_list) > 1):
            # 排除自交互和反向交互蛋白对
            if(PP_interaction_lines_list[0] == PP_interaction_lines_list[1] or
                    (PP_interaction_lines_list[1] + " " + PP_interaction_lines_list[0]) in PP_interaction):
                continue
            else:
                PP_interaction.append(PP_interaction_lines_list[0] + " " + PP_interaction_lines_list[1])

    for i in PP_interaction:
        PP_interaction_list = i.split(" ")
        if(PP_interaction_list[0] in all_p and PP_interaction_list[1] in all_p):
            # 无向连接图的邻接矩阵对称
            PM[all_p.index(PP_interaction_list[0])][all_p.index(PP_interaction_list[1])]=1
            PM[all_p.index(PP_interaction_list[1])][all_p.index(PP_interaction_list[0])] = 1

    # k=3 # Specifies the number of random walks
    random.choice(range(0, len(all_p)))


    # 获取随机游走过程中所有可能的noninteraction
    def get_all_possible_noninteraction(i):
        possible_noninteraction_protein=[]
        random_data1 = []
        for ind in range(len(all_p)):
            if (PM[i][ind] == 1):
                random_data1.append(ind)
        if (len(random_data1) > 0):
            for ind1 in random_data1:
                random_data2 = []
                for ind2 in range(len(all_p)):
                    if (PM[ind1][ind2] == 1):
                        random_data2.append(ind2)
                if(len(random_data2)>0):
                    for ind3 in random_data2:
                        random_data3 = []
                        for ind4 in range(len(all_p)):
                            if (PM[ind3][ind4] == 1):
                                random_data3.append(ind4)
                            if(len(random_data3)>0):
                                for ind_end in random_data3:
                                    if(all_p[ind_end] in possible_noninteraction_protein): pass
                                    elif(((all_p[i]+ " " +all_p[ind_end]) in PP_interaction) or
                                         ((all_p[ind_end]+ " " +all_p[i]) in PP_interaction)):
                                        pass
                                    else:possible_noninteraction_protein.append(all_p[ind_end])
        return possible_noninteraction_protein

    every_count=53
    for i in  range(len(all_p)):
        possible_noninteraction_protein_list=get_all_possible_noninteraction(i)
        if(len(possible_noninteraction_protein_list)>every_count):
            # 随机选择
            PPN=random.sample(possible_noninteraction_protein_list,every_count)
        else:
            PPN=possible_noninteraction_protein_list
        with open("PP_noninteraction.txt", "a") as f:
            for p in PPN:
                f.write(all_p[i])
                f.write(" ")
                f.write(p)
                f.write("\n")


if __name__ == '__main__':
    get_PP_noninteraction()




import numpy as np
import random


def prepare_ppi(all_ppi_file_path,all_protein_id_file_path):
    score=600  # Protein connectivity score
    save_list=[]
    with open(all_ppi_file_path, "r") as f:
        data_ppi = f.read()
        
    with open(all_protein_id_file_path, "r") as f:
        data_protein_id = f.read()
        
    data_ppi_list=data_ppi.split("\n")
    data_protein_id_list=data_protein_id.split("\n")
    for i in data_ppi_list:
        if(i!=""):
            if((i.split(" ")[0] in data_protein_id_list) and (i.split(" ")[1] in data_protein_id_list)
                    and (int(i.split(" ")[2])>=score) and (i.split(" ")[0]!=i.split(" ")[1])):
                temp = []
                temp.append(i.split(" ")[0])
                temp.append(i.split(" ")[1])
                save_list.append(temp)
                print("count + 1 . Already obtained " + str(len(save_list)))
    print("Total gain PPI interacting pair: " + str(len(save_list)))
    save_file_name="prepared_ppi_interacting_pair_moreThan"+str(score)+".txt"
    with open(save_file_name, "w") as f:
        for i in save_list:
            f.write(i[0])
            f.write(" ")
            f.write(i[1])
            f.write("\n")
    print("write success!!!, file name:"+save_file_name)


def ppi_txt_2_ppi_npy(ppi_file_path,all_protein_id_file_path):
    save_list=[]
    with open(ppi_file_path, "r") as f:
        data_ppi = f.read()
       
    with open(all_protein_id_file_path, "r") as f:
        data_protein_id = f.read()
        
    data_ppi_list = data_ppi.split("\n")
    data_protein_id_list = data_protein_id.split("\n")
    for i in data_ppi_list:
        if (i != ""):
            lines_list = i.split(" ")
            temp = []
            temp.append(data_protein_id_list.index(lines_list[0])+1)
            temp.append(data_protein_id_list.index(lines_list[1])+1)
            save_list.append(temp)

    save_path=ppi_file_path.split(".")[0]+".npy"
    np.save(save_path, save_list)



def prepare_protein_GO_interacting(protein_GO_file_path):
    save_list=[]
    with open(protein_GO_file_path, "r") as f:
        data_protein_GO = f.read()
        
    data_protein_GO_list=data_protein_GO.split("\n")
    for i in data_protein_GO_list:
        if(i!=""):
            i_list=i.split("\t")
            GO_list=i_list[1].split("; ")
            for j in GO_list:
                if(j!=""):
                    temp=[]
                    temp.append(i_list[0])
                    temp.append(j)
                    save_list.append(temp)
    with open("protein-GO_interacting_pair.txt", "w") as f:
        for i in save_list:
            f.write(i[0])
            f.write(" ")
            f.write(i[1])
            f.write("\n")
    print("write success!!!, file name: protein-GO_interacting_pair.txt")


# Get non interacting protein pairs
# all_ppi_file_path  is a file containing interactive protein pairs (this paper uses the PPI network of all human protein interactions downloaded from the string database),all_protein_id_file_path
# 输出
def get_non_interacting_protein(all_ppi_file_path,all_protein_id_file_path):
    max_every_pritein_interacting_count=53
    non_interacting_pair_txt=[]
    non_interacting_pair_npy=[]
    with open(all_ppi_file_path, "r") as f:
        data = f.read()
        
    with open(all_protein_id_file_path, "r") as f:
        data_protein_id = f.read()
        
    data_list=data.split("\n")
    data_protein_id_list=data_protein_id.split("\n")
    data_protein_id_list_random=data_protein_id.split("\n") #Reorder arrays randomly
    dict_protein_interacting={} # A dictionary of other proteins to which a protein is linked
    for i in data_protein_id_list:
        dict_protein_interacting[i]=[] #Initialize dictionary
    for i in data_list:
        if(i!=""):
            if i.split(" ")[0] in dict_protein_interacting.keys():
                dict_protein_interacting[i.split(" ")[0]].append(i.split(" ")[1]) #获取
            else:
                dict_protein_interacting[i.split(" ")[0]] = []
                dict_protein_interacting[i.split(" ")[0]].append(i.split(" ")[1])
    print("get protein interacting pair dictionary success！！")
    print("get non-interacting pair，start......")
    for i in data_protein_id_list:
        if(i!=""):
            print("getting "+i+" non-interacting pair")
            count = 0
            random.shuffle(data_protein_id_list_random)
            for j in data_protein_id_list_random:
                if((j in dict_protein_interacting[i]) or i==j or j==""):
                    pass
                else:
                    temp_txt=[i,j]
                    non_interacting_pair_txt.append(temp_txt)
                    temp_npy=[]
                    temp_npy.append(data_protein_id_list.index(i)+1)
                    temp_npy.append(data_protein_id_list.index(j)+1)
                    non_interacting_pair_npy.append(temp_npy)
                    count+=1
                    
                    if(count>=max_every_pritein_interacting_count):
                        print("non-interacting pair count == "+str(max_every_pritein_interacting_count)+"max count")
                        break
    
    print("Start writing。。。。non_interacting_pair_txt.txt")
    file_path_txt="non_interacting_pair_ppi_every"+str(max_every_pritein_interacting_count)
    with open(file_path_txt+".txt", "w") as f:
        for i in non_interacting_pair_txt:
            f.write(i[0])
            f.write(" ")
            f.write(i[1])
            f.write("\n")
    print("Start writing。。。non_interacting_pair_npy.npy")
    np.save(file_path_txt+".npy", non_interacting_pair_npy)


# According to all proteins_ Id.txt generates the corresponding protein sequence file, which will generate two files（1. Separated with all protein_ Id.txt corresponds to a file that contains only the protein sequence. TXT,, 2. An. NPY file of the corresponding protein sequence）
# Parameter file_ Path is the file path of. FA file of all human protein sequences exported from string database，all_protein_id_file_path
def save_all_seq(all_sequence_file_path,all_protein_id_file_path):
    pro_id_list=[]
    seq_list=[]
    save_seq_list=[]
    # npy_list=[]
    with open(all_sequence_file_path, "r") as f:
        data_sequence = f.read()
        
    with open(all_protein_id_file_path, "r") as f:
        data_protein_id = f.read()
        
    data_sequence_list=data_sequence.split(">")
    data_protein_id_list=data_protein_id.split("\n")


    for i in data_sequence_list:
        if(i!=""):
            i_list=i.split("\n")
            pro_id = i_list[0]
            del(i_list[0])
            pro_seq = ''.join(i_list)
            pro_id_list.append(pro_id)
            seq_list.append(pro_seq)

    for i in data_protein_id_list:
        if(i in pro_id_list):
            save_seq_list.append(seq_list[data_protein_id_list.index(i)])

    with open("protein_sequences_all.txt","w") as f:
        for i in save_seq_list:
            f.write(i)
            f.write("\n")
    np.save('protein_sequences_all.npy', save_seq_list)

def get_sequence_len(file_path):
    with open(file_path, "r") as f:
        data_sequence = f.read()
        
    data_sequence_list=data_sequence.split("\n")
    count=0
    more6000count=0
    for i in data_sequence_list:
        if(len(i)>count):
            count=len(i)
            print(count)
        if(len(i)>6000):
            more6000count+=1
    print(count)
    print("more than 6000")
    print(more6000count)

def graph2int(graph_file_path,node_file_path):
    node2int=[]
    with open(graph_file_path, "r") as f:
        data_graph = f.read()
        
    with open(node_file_path, "r") as f:
        data_node = f.read()
        
    data_graph_list=data_graph.split("\n")
    data_node_list=data_node.split("\n")
    for i in data_graph_list:
        if(i!=""):
            temp=[]
            temp.append(data_node_list.index(i.split("\t")[0])+1)
            temp.append(data_node_list.index(i.split("\t")[1])+1)
            node2int.append(temp)
    print("writing......")
    with open("graph2int_interacting_pair.txt","w") as f:
        for i in node2int:
            f.write(str(i[0]))
            f.write("\t")
            f.write(str(i[1]))
            f.write("\n")
    print("write success!!!")



def read_npy(file_path):
    
    arr=np.load(file_path)
    print(arr)
    # print(len(arr[0][0]))
    print(len(arr))
    print(type(arr))
    # print(type(arr[0]))
    # print(arr[0])
    # for i in arr[0]:
    #     print(i)

    # print(len(arr[0]))
    # print(len(arr[0][0]))
    print("load .npy done")


def print_hi(name):
    
    print(f'***** {name} *****')  


if __name__ == '__main__':
    print_hi('Hello world')
    # get_sequence_len("protein_sequences_all.txt")
    # prepare_ppi("9606.protein.links.v11.0.txt", "all-protein_id.txt")
    # prepare_protein_GO_interacting("protein-GO.txt")
    # get_non_interacting_protein("9606.protein.links.v11.0.txt", "all-protein_id.txt")
    # save_all_seq("9606.protein.sequences.v11.5.fa","all-protein_id.txt")
    # ppi_txt_2_ppi_npy("prepared_ppi_interacting_pair_.txt", "all-protein_id.txt")
    # read_npy("")
    # graph2int("proteinAndGO_all_interacting_pair.txt", "proteinAndGO_node.txt")

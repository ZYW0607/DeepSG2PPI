import numpy as np

def get_GO_features(ppi_file_path,node2vec_result_file_path,protein_seq_file_path):
    save_list=[]
    GO_features_list=[]

    data_ppi_list = np.load(ppi_file_path)
    sequence = np.load(protein_seq_file_path)
    with open(node2vec_result_file_path, "r") as f:
        data_node2vec_result = f.read()
        print("The file has " + str(len(data_node2vec_result)) + "lines in total")
        #print(type(data_node2vec_result))
   
    data_node2vec_result_list=data_node2vec_result.split("\n")
    # data_protein_seq_list=data_protein_seq.split("\n")

    for i in data_ppi_list:
        if(len(sequence[i[0]-1])<1800 and len(sequence[i[1]-1])<1800):
            temp=data_node2vec_result_list[i[0]-1].split("\t")
            temp.extend(data_node2vec_result_list[i[1]-1].split("\t"))
            numbers = list(map(float, temp)) # Converts a string in a list to a number
            save_list.append(numbers)
            print(len(save_list))

    print("Start writing...GO_features.npy")
    np.save("load_data_GO_features_interacting_.npy", save_list)
    # np.save("load_data_GO_features_non_interacting_.npy", save_list)


def read_npy(file_path):
    
    arr=np.load(file_path,allow_pickle=True)


    print("load .npy done")



if __name__ == '__main__':

    get_GO_features("ppi_interacting_pair_.npy","./node2vec_result/node2vec_result_only_vec.txt","protein_sequences_all.npy")
    # get_GO_features("non_interacting_pair_.npy","./node2vec_result/node2vec_result_only_vec.txt","protein_sequences_all.npy")
    # read_npy('GO_features600.npy')

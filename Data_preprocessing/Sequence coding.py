# import keras
'''
To generate 2-channel PPI data 
'''
import numpy as np

# A C D E F G H I K L M N P Q R S T V W Y



def get_CNN_data_2_channel(protein_A_seq,protein_B_seq):
    def trans_amio_acid_value(char):
        amio_acid_value = 0
        if char == 'A':
            amio_acid_value = 1.11
        elif char == 'C':
            amio_acid_value = 1.15
        elif char == 'D':
            amio_acid_value = 1.19
        elif char == 'E':
            amio_acid_value = 1.23
        elif char == 'F':
            amio_acid_value = 1.27
            # ---------------------------------
        elif char == 'G':
            amio_acid_value = 1.31
        elif char == 'H':
            amio_acid_value = 1.35
        elif char == 'I':
            amio_acid_value = 1.39
        elif char == 'K':
            amio_acid_value = 1.43
        elif char == 'L':
            amio_acid_value = 1.47
            # ---------------------------------
        elif char == 'M':
            amio_acid_value = 1.51
        elif char == 'N':
            amio_acid_value = 1.55
        elif char == 'P':
            amio_acid_value = 1.59
        elif char == 'Q':
            amio_acid_value = 1.63
        elif char == 'R':
            amio_acid_value = 1.67
            # ---------------------------------
        elif char == 'S':
            amio_acid_value = 1.71
        elif char == 'T':
            amio_acid_value = 1.75
        elif char == 'V':
            amio_acid_value = 1.79
        elif char == 'W':
            amio_acid_value = 1.83
        elif char == 'Y':
            amio_acid_value = 1.87
        else:
            amio_acid_value = 1.91

        return amio_acid_value
    def get_channel_1_sequence(seq_list):
        channel_data = []
        for i in range(len(seq_list)):
            channel_data.append(trans_amio_acid_value(seq_list[i]))
        assert len(seq_list) == len(channel_data),"error: wanted:{0} but got:{1}".format(len(seq_list),len(channel_data))
        return channel_data

    def get_channel_2_group(seq_list):

        value_list = []
        for i in range(len(seq_list)):
            if(i!=0 and i!=len(seq_list)-1):
                value_list.append(trans_amio_acid_value(seq_list[i-1])*2+(trans_amio_acid_value(seq_list[i+1]))*21)
            elif(i == 0):
                value_list.append(0*2+(trans_amio_acid_value(seq_list[i+1])) * 21)
            else:
                value_list.append(trans_amio_acid_value(seq_list[i-1])*2 + 0 * 21)


        # value_list.append(0)
        assert len(seq_list) == len(value_list), "error: wanted:{0} but got:{1}".format(len(seq_list),
                                                                                           len(value_list))

        return value_list

    channal_matrix = np.zeros((3600,2))

    c11 = get_channel_1_sequence(protein_A_seq)
    c13 = get_channel_2_group(protein_A_seq)

    c21 = get_channel_1_sequence(protein_B_seq)
    c23 = get_channel_2_group(protein_B_seq)

    for i in range(len(protein_A_seq)):
        channal_matrix[1799 - i][0] = c11[i]
        channal_matrix[1799 - i][1] = c13[i]
    for j in range(len(protein_B_seq)):
        channal_matrix[1800 + j][0] = c21[j]
        channal_matrix[1800 + j][1] = c23[j]
    return channal_matrix

def get_channal_matrix_features():


    map_list = np.load("interaction_pair_.npy")


    sequence = np.load("protein_sequences_all.npy")

    all_data = []
    for i in range(len(map_list)):

        s1 = sequence[int(map_list[i][0])].replace('\n','')

        s2 = sequence[int(map_list[i][1])].replace('\n','')

        if (len(s1)<1800) and ((len(s2)<1800)):
            all_data.append(get_CNN_data_2_channel(s1,s2))

    print("start save .npy")

    np.save("load_data_interacting.npy",all_data)

def get_protein_statistics_features(protein_A_seq,protein_B_seq):
    def get_protein_statistics(seq_list):
        total_data = [0 for i in range(21)]

        # A C D E F G H I K L M N P Q R S T V W Y ..
        for sta in seq_list:
            if sta == 'A':
                total_data[0] = total_data[0] + 1
            elif sta == 'C':
                total_data[1] = total_data[1] + 1
            elif sta == 'D':
                total_data[2] = total_data[2] + 1
            elif sta == 'E':
                total_data[3] = total_data[3] + 1
            elif sta == 'F':
                total_data[4] = total_data[4] + 1
            elif sta == 'G':
                total_data[5] = total_data[5] + 1
            elif sta == 'H':
                total_data[6] = total_data[6] + 1
            elif sta == 'I':
                total_data[7] = total_data[7] + 1
            elif sta == 'K':
                total_data[8] = total_data[8] + 1
            elif sta == 'L':
                total_data[9] = total_data[9] + 1
            elif sta == 'M':
                total_data[10] = total_data[10] + 1
            elif sta == 'N':
                total_data[11] = total_data[11] + 1
            elif sta == 'P':
                total_data[12] = total_data[12] + 1
            elif sta == 'Q':
                total_data[13] = total_data[13] + 1
            elif sta == 'R':
                total_data[14] = total_data[14] + 1
            elif sta == 'S':
                total_data[15] = total_data[15] + 1
            elif sta == 'T':
                total_data[16] = total_data[16] + 1
            elif sta == 'V':
                total_data[17] = total_data[17] + 1
            elif sta == 'W':
                total_data[18] = total_data[18] + 1
            elif sta == 'Y':
                total_data[19] = total_data[19] + 1
            else:
                total_data[20] = total_data[20] + 1

        for i in range(len(total_data)):
            total_data[i]=(total_data[i] / len(seq_list))

        return total_data

    protein_statistics_features = [0 for i in range(42)]
    s1=get_protein_statistics(protein_A_seq)
    s2=get_protein_statistics(protein_B_seq)
    for i in range(21):
        protein_statistics_features[20-i]=s1[i]
        protein_statistics_features[21+i]=s2[i]


    return protein_statistics_features

def get_all_features():

    map_list = np.load("interacting_pair.npy")


    sequence = np.load("protein_sequences_all.npy")
    
    data_channal_matrix = []
    data_statistics_features=[]
    for i in range(len(map_list)):

        s1 = sequence[int(map_list[i][0])-1].replace('\n', '')

        s2 = sequence[int(map_list[i][1])-1].replace('\n', '')

        if (len(s1) < 1800) and ((len(s2) < 1800)):
            data_channal_matrix.append(get_CNN_data_2_channel(s1, s2))
            # data_statistics_features.append(get_protein_statistics_features(s1, s2))

    print("start save seq_features .npy")


    np.save("load_data_interacting_matrix.npy", data_channal_matrix)



if __name__ == '__main__':
    get_all_features()
    # pass
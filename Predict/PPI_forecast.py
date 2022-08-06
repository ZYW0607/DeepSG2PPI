import numpy as np


#获取两通道的序列编码
def get_CNN_data_2_channel(protein_A_seq, protein_B_seq):

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
        assert len(seq_list) == len(channel_data), "error: wanted:{0} but got:{1}".format(len(seq_list),
                                                                                          len(channel_data))
        return channel_data

    def get_channel_2_group(seq_list):

        value_list = []
        for i in range(len(seq_list)):
            if (i != 0 and i != len(seq_list) - 1):
                value_list.append(
                    trans_amio_acid_value(seq_list[i - 1]) * 2 + (trans_amio_acid_value(seq_list[i + 1])) * 21)
            elif (i == 0):
                value_list.append(0 * 2 + (trans_amio_acid_value(seq_list[i + 1])) * 21)
            else:
                value_list.append(trans_amio_acid_value(seq_list[i - 1]) * 2 + 0 * 21)

        # value_list.append(0)
        assert len(seq_list) == len(value_list), "error: wanted:{0} but got:{1}".format(len(seq_list),
                                                                                        len(value_list))

        return value_list

    # if (len(protein_A_seq) < 1800) and ((len(protein_B_seq) < 1800)):
    sequence_features_list = []

    channal_matrix = np.zeros((3600, 2))

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



    sequence_features_list.append(channal_matrix)
    # print("得到 matrix features")
    return sequence_features_list



#获取统计特征
def get_protein_statistics_features(protein_A_seq, protein_B_seq):
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
            total_data[i] = (total_data[i] / len(seq_list))

        return total_data

    protein_statistics_features = [0 for i in range(42)]
    s1 = get_protein_statistics(protein_A_seq)
    s2 = get_protein_statistics(protein_B_seq)
    for i in range(21):
        protein_statistics_features[20 - i] = s1[i]
        protein_statistics_features[21 + i] = s2[i]

    protein_statistics_features_list=[]
    protein_statistics_features_list.append(protein_statistics_features)
    # print("得到 statistics features")
    return protein_statistics_features_list

#获取GO功能特征
def get_GO_features(protein_A_seq, protein_B_seq, node2vec_result_list, protein_seq_list):
    save_list = []

    # if (len(sequence[i[0] - 1]) < 1800 and len(sequence[i[1] - 1]) < 1800):
    # if (len(protein_A_seq) < 1800 and len(protein_B_seq) < 1800):
    temp = node2vec_result_list[protein_seq_list.index(protein_A_seq)].split("\t")
    temp.extend(node2vec_result_list[protein_seq_list.index(protein_B_seq)].split("\t"))
    numbers = list(map(float, temp))  # Converts a string in a list to a number
    save_list.append(numbers)
    # print(len(save_list))

    # print("得到GO_features.npy")
    # np.save("load_data_GO_features_interacting_.npy", save_list)
    # np.save("load_data_GO_features_non_interacting_.npy", save_list)
    return save_list


from keras.models import load_model

from numpy import argmax



from keras.utils.generic_utils import CustomObjectScope
import keras
from CA import _CA

import keras.backend as K
def matthews_correlation(y_true, y_pred):
    y_pred_pos = K.round(K.clip(y_pred, 0, 1))
    y_pred_neg = 1 - y_pred_pos

    y_pos = K.round(K.clip(y_true, 0, 1))
    y_neg = 1 - y_pos

    tp = K.sum(y_pos * y_pred_pos)
    tn = K.sum(y_neg * y_pred_neg)

    fp = K.sum(y_neg * y_pred_pos)
    fn = K.sum(y_pos * y_pred_neg)

    numerator = (tp * tn - fp * fn)
    denominator = K.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    return numerator / (denominator + K.epsilon())

def FPR(y_true, y_pred):
    y_pred_pos = K.round(K.clip(y_pred, 0, 1))
    y_pred_neg = 1 - y_pred_pos

    y_pos = K.round(K.clip(y_true, 0, 1))
    y_neg = 1 - y_pos

    tp = K.sum(y_pos * y_pred_pos)
    tn = K.sum(y_neg * y_pred_neg)

    fp = K.sum(y_neg * y_pred_pos)
    fn = K.sum(y_pos * y_pred_neg)

    return fp / (tn + fp)

def TPR(y_true, y_pred):
    y_pred_pos = K.round(K.clip(y_pred, 0, 1))
    y_pred_neg = 1 - y_pred_pos

    y_pos = K.round(K.clip(y_true, 0, 1))
    y_neg = 1 - y_pos

    tp = K.sum(y_pos * y_pred_pos)
    tn = K.sum(y_neg * y_pred_neg)

    fp = K.sum(y_neg * y_pred_pos)
    fn = K.sum(y_pos * y_pred_neg)

    return tp / (tp + fn)

def get_lr_metric(optimizer):  # printing the value of the learning rate
    def lr(y_true, y_pred):
        return optimizer.lr
    return lr
# lr = get_lr_metric(opt)
#
def get_lr():
    return 1
import keras_metrics as km
# import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import plot_model,to_categorical

import tensorflow as tf
config =tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
session = tf.compat.v1.InteractiveSession(config=config)



if __name__ == '__main__':
    with open("all-protein_id.txt", "r") as f:  # 打开文件
        protein_id = f.read()  # 读取文件
        # print(protein_id)

    with open("protein_sequences_all.txt", "r") as f:  # 打开文件
        protein_seq = f.read()  # 读取文件
        # print(protein_seq)

    protein_id_list = protein_id.split("\n")
    protein_seq_list = protein_seq.split("\n")

    # print(protein_id_list)
    # print(protein_id_list[0])
    # print(len(protein_id_list[0]))
    # print(protein_id_list.index(protein_id_list[0]))
    # print(protein_seq_list[protein_id_list.index(protein_id_list[0])])

    with open("node2vec_result_vec.txt", "r") as f:
        node2vec_result = f.read()
        print("The file has " + str(len(node2vec_result)) + "lines in total")
        # print(type(data_node2vec_result))

    node2vec_result_list = node2vec_result.split("\n")

    with open("Protein2predict.txt", "r") as f:  # 打开待预测蛋白质文件
        predicting = f.read()  # 读取文件
        print(protein_seq)
    predicting_list = predicting.split("\n")

    with CustomObjectScope(
            {'LeakyReLU': keras.layers.advanced_activations.LeakyReLU, 'TPR': TPR, 'FPR': FPR, 'lr': get_lr_metric,
             'binary_f1_score': km.f1_score(), 'matthews_correlation': matthews_correlation, 'ACCC': get_lr}):
        model = load_model('Trained_model.h5')  # 选取自己的.h5模型名称

    matrix_feature=""
    statistics_feature=""
    GO_feature=""
    # A=0
    if(len(predicting_list)>1):
        for A in range(0, len(predicting_list)):
            print("第"+str(A)+"个")
            # print(A)
            if(predicting_list[A] in protein_id_list):
                seqA =protein_seq_list[protein_id_list.index(predicting_list[A])]
                if (len(seqA) > 1800):
                    print("大于1800")
                if(len(seqA)< 1800):
                    for B in predicting_list:

                        if ((B in protein_id_list) and B!=predicting_list[A] and len(protein_seq_list[protein_id_list.index(B)])<1800):
                            seqB=protein_seq_list[protein_id_list.index(B)]
                            matrix_feature=get_CNN_data_2_channel(seqA,seqB)

                            matrix_feature = np.array(matrix_feature).reshape(len(matrix_feature), 60, 60, 2)
                            statistics_feature=get_protein_statistics_features(seqA,seqB)
                            statistics_feature = np.array(statistics_feature).reshape(len(statistics_feature), 42, 1)
                            GO_feature=get_GO_features(seqA,seqB,node2vec_result_list,protein_seq_list)
                            GO_feature = np.array(GO_feature).reshape(len(GO_feature), 256, 1)

                            pred_y = model.predict([matrix_feature, statistics_feature, GO_feature], verbose=0)
                            # print(pred_y)
                            predict = argmax(pred_y, axis=1)

                            if(predict[0]==1):
                                with open("PP_interaction.txt", "a") as f:
                                    f.write(protein_id_list[protein_seq_list.index(seqA)])
                                    f.write("\t")
                                    f.write(protein_id_list[protein_seq_list.index(seqB)])
                                    f.write("\n")
                            if (predict[0] == 0):
                                with open("PP_noninteraction.txt", "a") as f:
                                    f.write(protein_id_list[protein_seq_list.index(seqA)])
                                    f.write("\t")  
                                    f.write(protein_id_list[protein_seq_list.index(seqB)])
                                    f.write("\n")







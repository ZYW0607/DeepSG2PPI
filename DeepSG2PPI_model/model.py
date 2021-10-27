'''
 
'''

import os
import pickle
import keras
import numpy as np
from keras.models import Sequential
from keras.models import Model
from keras.layers import Dropout, Dense, Conv2D, Flatten, MaxPooling2D, Input, advanced_activations,Activation,Convolution1D,MaxPool1D,BatchNormalization
from keras.layers import Embedding, Conv1D, MaxPooling1D, GlobalMaxPooling1D,Reshape
from keras.callbacks import ModelCheckpoint
from keras import utils  
from sklearn.model_selection import train_test_split
from keras.layers.merge import concatenate
from keras.utils import plot_model,to_categorical
from keras.optimizers import Adam

import CBAM
import ECAnet
from CA import _CA


import tensorflow as tf
from sklearn.model_selection import train_test_split

import keras_metrics as km
import keras.backend as K
import pickle

config =tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
session = tf.compat.v1.InteractiveSession(config=config)




test_rate=0.1

def load_data_seq(file_path_interacting_seq,file_path_non_interacting_seq):
    interacting_data_seq= np.load(file_path_interacting_seq)
    non_interacting_data_seq = np.load(file_path_non_interacting_seq)

    y_train = [1 for s in range(len(interacting_data_seq)+len(non_interacting_data_seq))]
    print("y_train1111的len为：" + str(len(y_train)))
    y_train[len(interacting_data_seq):]=[0 for s in range(len(non_interacting_data_seq))]

    seq=np.concatenate((interacting_data_seq, non_interacting_data_seq), axis=0)

    random_seed=2
    train_x_seq, val_x_seq, train_y, val_y = train_test_split(seq, y_train, test_size=test_rate, random_state=random_seed,
                                                      shuffle=True)


    return train_x_seq, val_x_seq,train_y, val_y

def load_data_GO(file_path_interacting_GO,file_path_non_interacting_GO):

    interacting_data_GO = np.load(file_path_interacting_GO)
    non_interacting_data_GO = np.load(file_path_non_interacting_GO)


    GO=np.concatenate((interacting_data_GO, non_interacting_data_GO), axis=0)
    random_seed=2

    train_x_GO, val_x_GO = train_test_split(GO, test_size=test_rate, random_state=random_seed,
                                                              shuffle=True)


    return train_x_GO, val_x_GO

def load_data_seq_statistics(file_path_interacting_seq_statistics,file_path_non_interacting_seq_statistics):
    interacting_data_seq_statistics= np.load(file_path_interacting_seq_statistics)
    non_interacting_data_seq_statistics = np.load(file_path_non_interacting_seq_statistics)




    seq_statistics=np.concatenate((interacting_data_seq_statistics, non_interacting_data_seq_statistics), axis=0)

    random_seed=2
    train_x_seq_statistics, val_x_seq_statistics = train_test_split(seq_statistics, test_size=test_rate, random_state=random_seed,
                                                      shuffle=True)

    return train_x_seq_statistics, val_x_seq_statistics

# train_x_seq, val_x_seq,train_y, val_y = load_data_seq(
#     "load_data_channal_matrix_test.npy",
#     "load_data_channal_matrix_test.npy")
#
# train_x_GO, val_x_GO = load_data_GO(
#     "GO_features600_5.npy",
#     "GO_features600_5.npy")
#
# train_x_seq_statistics, val_x_seq_statistics = load_data_seq_statistics(
#     "load_data_statistics_test.npy",
#     "load_data_statistics_test.npy")



train_x_seq = np.array(train_x_seq).reshape(len(train_x_seq),60,60,2)
val_x_seq = np.array(val_x_seq).reshape(len(val_x_seq),60,60,2)

train_x_seq_statistics = np.array(train_x_seq_statistics).reshape(len(train_x_seq_statistics),42,1)
val_x_seq_statistics = np.array(val_x_seq_statistics).reshape(len(val_x_seq_statistics),42,1)

train_x_GO = np.array(train_x_GO).reshape(len(train_x_GO),256,1)
val_x_GO = np.array(val_x_GO).reshape(len(val_x_GO),256,1)



input_GO = Input(shape=train_x_GO.shape[1:])
input_seq = Input(shape=train_x_seq.shape[1:])
input_seq_statistics = Input(shape=train_x_seq_statistics.shape[1:])


leakrelu = advanced_activations.LeakyReLU(alpha=0.3)


def seq_model(input_data):
    # x = Dense(250, input_shape=input_data.shape[1:])(input_data)
    # x = Conv2D(32, (3, 3), padding='same')(input_data)
    # x = BatchNormalization()(x)
    # x = Activation(leakrelu)(x)
    #
    # x = Conv2D(32, (3, 3), padding='same')(x)
    # x = BatchNormalization()(x)
    # x = Activation(leakrelu)(x)
    # x = MaxPooling2D(pool_size=(2, 2))(x)
    x = Conv2D(64, (3, 3), padding='same')(input_data)
    x = BatchNormalization()(x)
    x = Activation(leakrelu)(x)
    x = Conv2D(64, (3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Activation(leakrelu)(x)
    x = MaxPooling2D(pool_size=(2, 2))(x)

    x = Conv2D(128, (3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Activation(leakrelu)(x)
    x = Conv2D(128, (3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Activation(leakrelu)(x)
    x = MaxPooling2D(pool_size=(2, 2))(x)
    # x = Dropout(0.25)(x)

    x = Conv2D(256, (3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Activation(leakrelu)(x)
    x = Conv2D(256, (3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Activation(leakrelu)(x)
    x = MaxPooling2D(pool_size=(2, 2))(x)
    # x = Conv2D(512, (3, 3), padding='same')(x)
    # x = BatchNormalization()(x)
    # x = Activation(leakrelu)(x)
    # x = Conv2D(512, (3, 3), padding='same')(x)
    # x = BatchNormalization()(x)
    # x = Activation(leakrelu)(x)
    # x = MaxPooling2D(pool_size=(2, 2))(x)
    # x = Dropout(0.25)(x)
    # x = CBAM.cbam_module(x)
    # x = ECAnet.eca_layer(x, num=1)
    x = _CA(x, "CA")
    x = Flatten()(x)
    seq_model_ = Model(inputs=input_data, outputs=x)
    return seq_model_


def seq_statistics_model(input_seq_statistics):
    seq_statistics_model_ = Sequential()
    seq_statistics_model_.add(Convolution1D(filters=32, kernel_size=3, strides=1, padding='same', input_shape=(42, 1)))
    seq_statistics_model_.add(BatchNormalization())
    seq_statistics_model_.add(Activation(leakrelu))
    seq_statistics_model_.add(Convolution1D(filters=32, kernel_size=3, strides=1, padding='same'))
    seq_statistics_model_.add(BatchNormalization())
    seq_statistics_model_.add(Activation(leakrelu))
    seq_statistics_model_.add(MaxPool1D(pool_size=5, strides=1, padding="valid"))
    seq_statistics_model_.add(Convolution1D(filters=16, kernel_size=3, strides=1, padding='same'))
    seq_statistics_model_.add(BatchNormalization())
    seq_statistics_model_.add(Activation(leakrelu))
    seq_statistics_model_.add(Convolution1D(filters=16, kernel_size=3, strides=1, padding='same'))
    seq_statistics_model_.add(BatchNormalization())
    seq_statistics_model_.add(Activation(leakrelu))
    
    seq_statistics_model_.add(MaxPool1D(pool_size=5, strides=1, padding="valid"))

    seq_statistics_model_.add(Flatten())

    return seq_statistics_model_


def GO_model(input_GO):
    GO_model_ = Sequential()
    GO_model_.add(Convolution1D(filters=64, kernel_size=5, strides=1, padding='same', input_shape=(256, 1)))
    GO_model_.add(BatchNormalization())
    GO_model_.add(Activation(leakrelu))
    GO_model_.add(Convolution1D(filters=64, kernel_size=5, strides=1, padding='same'))
    GO_model_.add(BatchNormalization())
    GO_model_.add(Activation(leakrelu))
    GO_model_.add(MaxPool1D(pool_size=5, strides=1, padding="valid"))
    GO_model_.add(Convolution1D(filters=32, kernel_size=5, strides=1, padding='same'))
    GO_model_.add(BatchNormalization())
    GO_model_.add(Activation(leakrelu))
    GO_model_.add(Convolution1D(filters=32, kernel_size=5, strides=1, padding='same'))
    GO_model_.add(BatchNormalization())
    GO_model_.add(Activation(leakrelu))
    GO_model_.add(MaxPool1D(pool_size=5, strides=1, padding="valid"))
    GO_model_.add(Flatten())

    return GO_model_


seq_model = seq_model(input_seq)
seq_statistics_model=seq_statistics_model(input_seq_statistics)
GO_model = GO_model(input_GO)

class_num=2
c = concatenate([seq_model.output, seq_statistics_model.output, GO_model.output], axis=-1)
z = Dense(1024)(c)
z = Activation(leakrelu)(z)
z = Dropout(0.3)(z)
z = Dense(1024)(z)
z = Activation(leakrelu)(z)
z = Dropout(0.3)(z)
z = Dense(class_num)(z)
z = Activation('softmax')(z)


model = Model(inputs=[seq_model.input, seq_statistics_model.input, GO_model.input], outputs=z)


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



opt = Adam(lr=0.001, decay=1e-3 / 100)


def get_lr_metric(optimizer):  # printing the value of the learning rate
    def lr(y_true, y_pred):
        return optimizer.lr
    return lr
lr = get_lr_metric(opt)
#
# from keras.callbacks import LearningRateScheduler
# def scheduler(epoch):
#     # Every 10 epochs, the learning rate is reduced to 1 / 10 of the original
#     if epoch % 10 == 0 and epoch != 0:
#         lr = K.get_value(model.optimizer.lr)
#         K.set_value(model.optimizer.lr, lr * 0.2)
#         print("lr changed to {}".format(lr * 0.2))
#     return K.get_value(model.optimizer.lr)
#
#
# reduce_lr = LearningRateScheduler(scheduler)

model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['binary_accuracy',lr,'mae','AUC','Precision','Recall',TPR,FPR,km.f1_score(),matthews_correlation])

model.summary()

"""Draw model drawing"""
plot_model(model, 'attention_&-1&input3_model_cbam.bmp', show_shapes=True)

train_y = to_categorical(train_y,2)
val_y = to_categorical(val_y,2)


np.random.seed(100)
np.random.shuffle(train_x_seq)
np.random.seed(100)
np.random.shuffle(train_x_GO)
np.random.seed(100)
np.random.shuffle(train_x_seq_statistics)
np.random.seed(100)
np.random.shuffle(train_y)


history=model.fit([train_x_seq, train_x_seq_statistics, train_x_GO], train_y,
          # validation_split=0.2,
          #validation_data=([val_x_seq, val_x_seq_statistics, val_x_GO], val_y),
          epochs=100, batch_size=60)
                  # ,callbacks=[reduce_lr])


# scores = model.evaluate([val_x_seq, val_x_seq_statistics, val_x_GO], val_y, verbose=0)
# 
# metrics_names=model.metrics_names
# evaluate_result=[]
# evaluate_result.append(metrics_names)
# evaluate_result.append(scores)
# file1 = open('./model_history/evaluate_result_filters_reduce.pkl', 'wb')
# pickle.dump(evaluate_result, file1)
# file1.close()

# 
pred_y = model.predict([val_x_seq, val_x_seq_statistics, val_x_GO], verbose=0)
y_val_pred={"pred_y":pred_y,"val_y":val_y}
file1 = open('./model_history/val_y&pred_y.pkl', 'wb')
pickle.dump(y_val_pred, file1)
file1.close()

# 

file = open('./model_history/history.pkl', 'wb')
pickle.dump(history.history, file)
file.close()

print("starting save model")
model.save('./model_history/model_.h5')



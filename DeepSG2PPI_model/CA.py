from tensorflow.python.keras.layers import Lambda, Concatenate, Reshape, Conv2D, BatchNormalization, Activation, \
    Multiply
import tensorflow as tf


def _CA(inputs, name, ratio=8):
    w, h, out_dim = [int(x) for x in inputs.shape[1:]]
    temp_dim = max(int(out_dim // ratio), ratio)

    h_pool = Lambda(lambda x: tf.reduce_mean(x, axis=1))(inputs)
    w_pool = Lambda(lambda x: tf.reduce_mean(x, axis=2))(inputs)

    x = Concatenate(axis=1)([h_pool, w_pool])
    x = Reshape((1, w + h, out_dim), name=name + '_Reshape')(x)
    x = Conv2D(temp_dim, 1)(x)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x_h, x_w = Lambda(lambda x: tf.split(x, [h, w], axis=2))(x)
    x_w = Reshape((w, 1, temp_dim))(x_w)

    x_w = Conv2D(out_dim, 1, activation='sigmoid')(x_w)
    x_h = Conv2D(out_dim, 1, activation='sigmoid')(x_h)
    x = Multiply()([inputs, x_h, x_w])
    return x
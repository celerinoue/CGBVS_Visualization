import tensorflow as tf
import numpy as np
import joblib
import kgcn.layers
import tensorflow.keras.layers as klayer
import tensorflow.contrib.keras as K
from kgcn.default_model import DefaultModel

class MultimodalNetwork(DefaultModel):
    def build_placeholders(self,info,config,batch_size,**kwargs):
        # input data types (placeholders) of this neural network
        return self.get_placeholders(info,config,batch_size,
            ['adjs','nodes','labels','mask','dropout_rate',
            'enabled_node_nums','is_train','features',
            'sequences','sequences_len','embedded_layer'],**kwargs)

    def build_model(self,placeholders,info,config,batch_size,feed_embedded_layer=False,**kwargs):
        self.batch_size = batch_size
        self.input_dim = info.input_dim
        self.adj_channel_num = info.adj_channel_num
        self.sequence_symbol_num = info.sequence_symbol_num
        self.graph_node_num = info.graph_node_num
        self.label_dim = info.label_dim
        self.embedding_dim = config["embedding_dim"]
        self.pos_weight = info.pos_weight
        self.feed_embedded_layer = feed_embedded_layer
        ## aliases
        batch_size = self.batch_size
        in_adjs = self.placeholders["adjs"]
        features = self.placeholders["features"]
        sequences = self.placeholders["sequences"]
        sequences_len = self.placeholders["sequences_len"]
        in_nodes = self.placeholders["nodes"]
        labels = self.placeholders["labels"]
        mask = self.placeholders["mask"]
        is_train = self.placeholders["is_train"]
        enabled_node_nums = self.placeholders["enabled_node_nums"]
        embedded_layer = self.placeholders['embedded_layer']

        layer = features

        #=== GCN =====
        print("graph input layer:",layer.shape)

        out_dim = 64
        layer = kgcn.layers.GraphConv(out_dim, self.adj_channel_num)(layer, adj=in_adjs)
        layer = kgcn.layers.GraphBatchNormalization()(layer, max_node_num=self.graph_node_num, enabled_node_nums=enabled_node_nums)
        layer = kgcn.layers.GraphDense(out_dim)(layer)
        layer = tf.nn.relu(layer)
        ###
        out_dim = 64
        layer = kgcn.layers.GraphConv(out_dim, self.adj_channel_num)(layer, adj=in_adjs)
        layer = kgcn.layers.GraphBatchNormalization()(layer, max_node_num=self.graph_node_num, enabled_node_nums=enabled_node_nums)
        layer = kgcn.layers.GraphDense(out_dim)(layer)
        layer = tf.nn.relu(layer)
        ###
        out_dim = 64
        layer = kgcn.layers.GraphConv(out_dim, self.adj_channel_num)(layer, adj=in_adjs)
        layer = kgcn.layers.GraphBatchNormalization()(layer, max_node_num=self.graph_node_num, enabled_node_nums=enabled_node_nums)
        layer = kgcn.layers.GraphDense(out_dim)(layer)
        layer = tf.nn.relu(layer)
        ###
        out_dim = 32
        layer = kgcn.layers.GraphDense(out_dim)(layer)
        layer = tf.nn.relu(layer)
        layer = kgcn.layers.GraphGather()(layer)
        graph_output_layer = layer

        print("graph output layer:",graph_output_layer.shape)

        #=== CNN =====

        with tf.variable_scope("seq_nn") as scope_part:
            # Embedding
            self.embedding_layer = tf.keras.layers.Embedding(self.sequence_symbol_num,self.embedding_dim)(sequences)

            if self.feed_embedded_layer:
                layer = embedded_layer
            else:
                layer = self.embedding_layer

            print("sequence input layer:",layer.shape)

            # CNN + Pooling
            # layer1
            layer = tf.keras.layers.Conv1D(filters=100, kernel_size=5, strides=1, padding="same", activation='relu')(layer)
            layer = tf.keras.layers.MaxPooling1D(pool_size=4)(layer)

            # layer2
            layer = tf.keras.layers.Conv1D(filters=100, kernel_size=5, strides=1, padding="same", activation='relu')(layer)
            layer = tf.keras.layers.MaxPooling1D(pool_size=4)(layer)

            # layer3
            layer = tf.keras.layers.Conv1D(filters=100, kernel_size=5, strides=1, padding="same", activation='relu')(layer)
            layer = tf.keras.layers.MaxPooling1D(pool_size=2)(layer)

            # down dimension
            layer = tf.keras.layers.Conv1D(filters=1, kernel_size=1, strides=1, padding="same", activation='tanh')(layer)
            layer = tf.squeeze(layer)
            if len(layer.shape) == 1:
                # When batch-size is 1, this shape doesn't have batch-size dimmension due to just previous 'tf.squeeze()'.
                layer = tf.expand_dims(layer, axis=0)
            seq_output_layer = layer

            print("sequence output layer:",seq_output_layer.shape)

        #=== shred part =====

        layer = tf.concat([seq_output_layer, graph_output_layer], axis=1)
        print("shared_part input:",layer.shape)

        with tf.variable_scope("shared_nn") as scope_part:
            layer = tf.keras.layers.BatchNormalization()(layer)
            layer = tf.keras.layers.Dense(64)(layer)
            layer = tf.keras.layers.BatchNormalization()(layer)
            layer = tf.nn.relu(layer)

        layer = tf.keras.layers.Dense(self.label_dim)(layer)

        # compute prediction
        prediction=tf.nn.softmax(layer)

        # compute cost and metrics
        cost=mask*tf.nn.softmax_cross_entropy_with_logits(labels=labels,logits=layer)
        cost_opt=tf.reduce_mean(cost)

        metrics={}
        cost_sum=tf.reduce_sum(cost)

        correct_count=mask*tf.cast(tf.equal(tf.argmax(prediction,1), tf.argmax(labels,1)),tf.float32)
        metrics["correct_count"]=tf.reduce_sum(correct_count)
        self.out=layer
        return self, prediction, cost_opt, cost_sum, metrics

    def embedding(self, sess, data):
        key = self.placeholders['sequences']
        feed_dict = {key: data}
        out = sess.run(self.embedding_layer, feed_dict)
        return out

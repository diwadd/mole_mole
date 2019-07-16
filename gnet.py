import tensorflow as tf


def weight_variable(shape, name):

    initial = tf.truncated_normal(shape, stddev=0.1)
    # tf.summary.histogram(nm, initial, collections=['always'])
    return tf.Variable(initial, name=name)


v1 = weight_variable([3, 3], "v1")


class GNet:

    def __init__(self, input_dim, state_dim, output_dim):

        self.input_dim = input_dim
        self.state_dim = state_dim
        self.output_dim = output_dim
        self.state_input = self.input_dim - 1 + state_dim

        self.state_l1 = 15
        self.state_l2 = self.state_dim
        self.output_l1 = 10
        self.output_l2 = self.output_dim

        self.weights = {'State_L1': weight_variable([self.state_input, self.state_l1], "WEIGHT_STATE_L1"),
                        'State_L2': weight_variable([ self.state_l1, self.state_l2], "WEIGHT_STATE_L1"),

                        'Output_L1':weight_variable([self.state_l2,self.output_l1], "WEIGHT_OUTPUT_L1"),
                        'Output_L2': weight_variable([self.output_l1, self.output_l2], "WEIGHT_OUTPUT_L2")
                        }

        self.biases = {'State_L1': weight_variable([self.state_l1],"BIAS_STATE_L1"),
                       'State_L2': weight_variable([self.state_l2], "BIAS_STATE_L2"),

                       'Output_L1':weight_variable([self.output_l1],"BIAS_OUTPUT_L1"),
                       'Output_L2': weight_variable([ self.output_l2], "BIAS_OUTPUT_L2")
                       }


net = GNet(3, 3, 3)


with tf.Session() as sess:
    init_op = tf.initialize_all_variables()
    sess.run(init_op)

    x = sess.run(net.weights)

    for k, v in x.items():
        print(k)
        print(v.shape)



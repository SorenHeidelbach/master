import numpy as np
from tensorflow import keras as ke
from tensorflow.keras import layers
import pandas as pd
import json

class AAE():
    def __init__(self, input_shape=(500,), latent_dim=2, hidden_nodes=1000, drop_out_rate=0.2, seed=2):
        self.latent_dim = latent_dim
        self.hidden_nodes = hidden_nodes
        self.drop_out_rate = drop_out_rate
        self.input_shape = input_shape
        self.initializer = ke.initializers.RandomNormal(mean=0.5, stddev=0.1, seed=seed)
        self.initialise()

    def _build_model_enc(self):
        input = ke.Input(shape=self.input_shape)
        x = layers.Dense(self.hidden_nodes, activation="relu", kernel_initializer=self.initializer)(input)
        x = layers.Dropout(self.drop_out_rate)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dense(self.hidden_nodes, activation="relu", kernel_initializer=self.initializer)(x)
        x = layers.BatchNormalization()(x)
        output = layers.Dense(self.latent_dim, name="encoder_out")(x)

        model = ke.Model(inputs=input, outputs=output, name="encoder")
        return model

    def _build_model_dec(self):
        input = ke.Input(shape=(self.latent_dim,), name="decoder_in")
        x = layers.Dense(self.hidden_nodes, activation="relu", kernel_initializer=self.initializer)(input)
        x = layers.BatchNormalization()(x)
        x = layers.Dense(self.hidden_nodes, activation="relu", kernel_initializer=self.initializer)(x)
        x = layers.BatchNormalization()(x)
        output = layers.Dense(self.input_shape[0])(x)

        model = ke.Model(inputs=input, outputs=output, name="decoder")
        return model

    def _build_model_disc(self):
        input = ke.Input(shape=(self.latent_dim,))
        x = layers.Dense(self.hidden_nodes, activation="relu", kernel_initializer=self.initializer)(input)
        x = layers.Dense(self.hidden_nodes, activation="relu", kernel_initializer=self.initializer)(x)
        output = layers.Dense(1, activation="sigmoid")(x)

        model = ke.Model(inputs=input, outputs=output, name="discriminator")
        return model

    def initialise(self):
        self.model_enc = self._build_model_enc()
        self.model_dec = self._build_model_dec()
        self.model_disc = self._build_model_disc()

        # Autoencoder
        input = ke.Input(shape=self.input_shape)
        encoded = self.model_enc(input)
        decoded = self.model_dec(encoded)
        self.model_ae = ke.Model(input, decoded, name="autoencoder")

        # Discriminator
        discrim = self.model_disc(encoded)
        self.model_enc_disc = ke.Model(input, discrim, name="encoder_discriminator")

        self.model_disc.compile(optimizer=ke.optimizers.Adam(learning_rate=1e-3), loss="binary_crossentropy")
        self.model_enc_disc.compile(optimizer=ke.optimizers.Adam(learning_rate=1e-3), loss="binary_crossentropy")
        self.model_ae.compile(optimizer=ke.optimizers.Adam(learning_rate=1e-3), loss="mse")

    def settrainable(self, model, toset):
        for layer in model.layers:
            layer.trainable = toset
        model.trainable = toset

    def train(self, x_train, batchsize=100, epochs=100):
        with open('/shared-nfs/SH/analysis_output/MGM1/training_log.txt', 'w') as f:
            f.write('reconstruction_loss\tadvesarial_loss\n')
        for epochnumber in range(epochs):
            np.random.shuffle(x_train)
            self.settrainable(self.model_ae, True)
            self.settrainable(self.model_enc, True)
            self.settrainable(self.model_dec, True)
            self.settrainable(self.model_enc_disc, True)
            for i in range(int(len(x_train) / batchsize)):
                batch = x_train[i * batchsize:(i + 1) * batchsize, :]
                self.model_ae.train_on_batch(batch, batch)

                self.settrainable(self.model_disc, True)
                batchpred = self.model_enc.predict(batch)
                fakepred = np.random.multivariate_normal([0] * self.latent_dim, np.identity(self.latent_dim), batchsize)
                discbatch_x = np.concatenate([batchpred, fakepred])
                discbatch_y = np.concatenate([np.zeros(batchsize), np.ones(batchsize)])
                self.model_disc.train_on_batch(discbatch_x, discbatch_y)

                self.settrainable(self.model_disc, False)
                self.model_enc_disc.train_on_batch(batch, np.ones(batchsize))
            print("Epoch number:", epochnumber)
            with open('/shared-nfs/SH/analysis_output/MGM1/training_log.txt', 'a') as f:
                f.write('{0}\t{1}\t{2}\n'.format(epochnumber,self.model_ae.evaluate(x_train, x_train, verbose=0), self.model_enc_disc.evaluate(x_train, np.ones(len(x_train)), verbose=0)))
            if epochnumber % 10 == 0:
                print("Reconstruction Loss:", self.model_ae.evaluate(x_train, x_train, verbose=0), "Adverserial Loss:", self.model_enc_disc.evaluate(x_train, np.ones(len(x_train)), verbose=0))
            if epochnumber % 100 == 0:
                np.savetxt("/shared-nfs/SH/analysis_output/MGM1/aae_epoch" + str(epochnumber) + ".tsv", self.model_enc.predict(x_train), delimiter="\t")


if __name__ == '__main__':

    train = pd.read_csv("/shared-nfs/SH/analysis_output/MGM1/features.csv").values
    x_train = train[:, 1:]
    x_train = x_train.astype(float)

    aae = AAE(input_shape=(x_train.shape[1],), hidden_nodes=100, latent_dim=10)
    aae.train(x_train, epochs=1000, batchsize=200)
    encoding = aae.model_enc.predict(x_train)
    np.savetxt("/shared-nfs/SH/analysis_output/MGM1/aae.tsv", encoding, delimiter="\t")

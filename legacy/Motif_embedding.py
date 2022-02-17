import numpy as np
from tensorflow import keras as ke
from tensorflow.keras import layers
import pandas as pd

class convAE():
    def __init__(self, input_shape=(7, 1), kernel_size = 3, stride = 1, latent_dim=2, drop_out_rate=0.2, seed=2):
        self.latent_dim = latent_dim
        self.kernel_size = kernel_size
        self.stride = stride
        self.drop_out_rate = drop_out_rate
        self.input_shape = input_shape
        self.initializer = ke.initializers.RandomNormal(mean=0.5, stddev=0.1, seed=seed)
        self.n_filter_1 = np.ceil((self.input_shape[0] - self.kernel_size + 1) / self.stride)
        self.n_filter_2 = np.ceil((self.n_filter_1 - self.kernel_size + 1) / self.stride)
        self.n_filter_3 = np.ceil((self.n_filter_2 - self.kernel_size + 1) / self.stride)
        self.initialise()

    def _build_model_enc(self):
        input = ke.Input(shape=self.input_shape)
        x = layers.Conv1D(self.n_filter_1, kernel_size=self.kernel_size)(input)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
        x = layers.Conv1D(self.n_filter_2, kernel_size=self.kernel_size)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
        x = layers.Conv1D(self.n_filter_3, kernel_size=2, activation="relu")(x)
        x = layers.Flatten()(x)
        output = layers.Dense(self.latent_dim, name="encoder_out")(x)

        model = ke.Model(inputs=input, outputs=output, name="encoder")
        return model

    def _build_model_dec(self):
        input = ke.Input(shape=(self.latent_dim,))
        x = layers.Reshape(target_shape=(2, 1))(input)
        x = layers.Conv1DTranspose(1, kernel_size=2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
        x = layers.Conv1DTranspose(1, kernel_size=self.kernel_size)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
        output = layers.Conv1DTranspose(1, kernel_size=self.kernel_size)(x)

        model = ke.Model(inputs=input, outputs=output, name="decoder")
        return model

    def initialise(self):
        self.model_enc = self._build_model_enc()
        self.model_dec = self._build_model_dec()

        # Autoencoder
        input = ke.Input(shape=self.input_shape)
        encoded = self.model_enc(input)
        decoded = self.model_dec(encoded)
        self.model_ae = ke.Model(input, decoded, name="autoencoder")
        self.model_ae.compile(optimizer=ke.optimizers.Adam(learning_rate=1e-3), loss="binary_crossentropy")

if __name__ == "__main__":
    train = pd.read_csv("/user/student.aau.dk/sheide17/motif_features.csv").values
    x_train = train.astype(float)
    x_train = x_train[:, :, np.newaxis]
    convae = convAE()
    convae.model_ae.summary()
    convae.model_ae.fit(
        x=x_train,
        y=x_train,
        epochs=20,
        batch_size=20000,
        shuffle=True,
        validation_data=(x_train, x_train),
    )
    endoded_output = convae.model_enc(x_train)
    np.savetxt("/shared-nfs/SH/analysis_output/MGM1/motif_convae.tsv", endoded_output.numpy(), delimiter="\t")


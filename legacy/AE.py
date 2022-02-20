import numpy as np
from tensorflow import keras as ke
from tensorflow.keras import layers
import pandas as pd
import os
import time
import json
from datetime import datetime


class AE():
    def __init__(self, input_shape=(500,), latent_dim=32, hidden_nodes=1000, learning_rate=1e-2, hidden_activation="relu"):
        # Encoder
        ## input
        enc_in = ke.Input(shape=input_shape)
        ## Hidden layer 1
        enc_x = layers.Dense(hidden_nodes, activation=hidden_activation)(enc_in)
        enc_x = layers.BatchNormalization()(enc_x)
        ## Hidden layer 2
        enc_x = layers.Dense(hidden_nodes, activation=hidden_activation)(enc_x)
        enc_x = layers.BatchNormalization()(enc_x)
        ## Latent representation
        enc_out = layers.Dense(latent_dim, activation="tanh")(enc_x)
        self.model_enc = ke.Model(inputs=enc_in, outputs=enc_out, name="encoder")

        # Decoder
        ## input, same dimension as latent representation
        dec_in = ke.Input(shape=(latent_dim,))
        ## Hidden layer 3
        dec_x = layers.Dense(hidden_nodes, activation=hidden_activation)(dec_in)
        dec_x = layers.BatchNormalization()(dec_x)
        ## Hidden layer 4
        dec_x = layers.Dense(hidden_nodes, activation=hidden_activation)(dec_x)
        dec_x = layers.BatchNormalization()(dec_x)
        ## Output layer
        dec_out = layers.Dense(input_shape[0])(dec_x)
        self.model_dec = ke.Model(inputs=dec_in, outputs=dec_out, name="decoder")

        # Autoencoder
        input = ke.Input(shape=input_shape)
        self.model_ae = ke.Model(input, self.model_dec(self.model_enc(input)), name="autoencoder")

        # Compile model
        self.model_ae.compile(optimizer=ke.optimizers.Adam(learning_rate=learning_rate), loss="mse")


if __name__ == "__main__":
    # Load settings
    with open("/shared-nfs/SH/code/params.json") as f:
        p = json.load(f)
    # Create folder for training dump
    save_dir = os.path.join(p["misc"]["output_folder"],
                             "AE_train_binning_" + datetime.fromtimestamp(time.time()).strftime("%d-%m-%y"))
    try:
        os.makedirs(save_dir)
    except FileExistsError:
        print("Directionary already exist for this training run, results will be overwritten")
        pass

    # Save settings for training run
    with open(save_dir + "/ae_settings.json", 'w') as f_json:
        json.dump(p["autoencoder"], f_json, indent = 4, sort_keys=True)

    # Load contig features
    features = pd.read_csv(p["misc"]["selected_features"]).values
    features = features[:, 1:] # drop contig names
    features = features.astype("float")

    # Initialise model
    ae = AE(
        input_shape=(features.shape[1],),
        hidden_nodes=p["autoencoder"]["intermediate_nodes"],
        latent_dim=p["autoencoder"]["bottleneck_nodes"],
        hidden_activation=p["autoencoder"]["activation_function"],
        learning_rate=p["autoencoder"]["learning_rate"]
    )

    # Train model
    ae.model_ae.fit(
        x=features,
        y=features,
        epochs=p["autoencoder"]["epochs"]
    )

    # Latent representation
    encoding = ae.model_enc.predict(features)

    # Save results
    np.savetxt(save_dir +"/binnning_AE_representation.tsv", encoding, delimiter="\t")
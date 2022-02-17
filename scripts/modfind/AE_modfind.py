#!/usr/bin/env python3

import argparse
# Parse arguments
# Instantiate the parser
parser = argparse.ArgumentParser(description='Embed multidimentional faeture matrix into fewer latent representations with autoencoder')
# Required positional arguments
parser.add_argument('feature_file', metavar="features",
                    help='TSV file containing abnomality feature matrix')
parser.add_argument('out',
                    help='Directionary to save results too')
# Optional arguments
parser.add_argument('--latent_nodes', type=int, default=6, metavar='latent',
                    help="Number of latent nodes in bottleneck layer")
parser.add_argument('--hidden_nodes', type=int, default=10, metavar='hidden',
                    help="Number of nodes in the hidden layers of the encoder and decoder")
parser.add_argument('--kernel_size', type=int, default=2, metavar='kernel_size',
                    help="Size of kernel in convolution layer")
parser.add_argument('-t', '--threads', type=int, default=10,
                    help="Threads to allocate training of autoencoder")
parser.add_argument('--epochs', type=int, default=10,
                    help="Epochs to train autoencoder for")
parser.add_argument('--rate', type=int, default=0.001, metavar='rate',
                    help="Training rate")
parser.add_argument('--batch_size', type=int, default=100, metavar='batch',
                    help="Size of the batches in each epoch")
parser.add_argument('--seed', type=int, default=42,
                    help="Seed for pytorch and numpy")
parser.add_argument('--val_percent', type=float, default=0.05,
                    help="Percent of data used for validation during training")
parser.add_argument('--event_frame_size', type=int, default=10,
                    help="Size of event frame used in preprocessing")

# Load arguments
arg = parser.parse_args()

import torch
from torch import nn
#torch.manual_seed(arg.seed)
#torch.set_num_threads(arg.threads)
import pandas as pd
import numpy as np
#np.random.seed(arg.seed)
import numpy as np
import pandas as pd
import os, json, csv, math


class Reshape(nn.Module):
    def __init__(self, *args):
        super().__init__()
        self.shape = args

    def forward(self, x):
        return x.view(self.shape)

class AE_conv(nn.Module):
    def __init__(self, input_shape, latent_nodes, hidden_nodes, kernel, stride):
        super().__init__()
        self.size_conv=input_shape[2] - 3 * (kernel - 1)
        self.encoder = nn.Sequential(
            nn.Conv1d(in_channels=2, out_channels=4, kernel_size=kernel, stride=stride),
            nn.BatchNorm1d(4),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            nn.Conv1d(in_channels=4, out_channels=8, kernel_size=kernel, stride=stride),
            nn.BatchNorm1d(8),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            nn.Conv1d(in_channels=8, out_channels=16, kernel_size=kernel, stride=stride),
            nn.BatchNorm1d(16),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            nn.Flatten(),
            nn.Linear(self.size_conv*16, latent_nodes),
            nn.BatchNorm1d(latent_nodes),
            nn.LeakyReLU()
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_nodes, hidden_nodes),
            nn.BatchNorm1d(hidden_nodes),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            nn.Linear(hidden_nodes, self.size_conv*16),
            nn.BatchNorm1d(self.size_conv*16),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            Reshape(-1, 16, self.size_conv),
            nn.ConvTranspose1d(in_channels=16, out_channels=8, kernel_size=kernel, stride=stride),
            nn.BatchNorm1d(8),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            nn.ConvTranspose1d(in_channels=8, out_channels=4, kernel_size=kernel, stride=stride),
            nn.BatchNorm1d(4),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.2),
            nn.ConvTranspose1d(in_channels=4, out_channels=2, kernel_size=kernel, stride=stride)
        )

    def forward(self, x, get_encoded=False):
        if get_encoded:
            x = self.encoder(x)
        else:
            x = self.encoder(x)
            x = self.decoder(x)
        return x

    def train(self, loader_train, loader_validate, epochs, learning_rate):
        # Validation using MSE Loss
        loss_function = torch.nn.MSELoss()

        # Initialise optimser
        optimiser = torch.optim.Adam(
            self.parameters(),
            lr=learning_rate,
            weight_decay=1e-8
        )

        outputs = []
        losses = []
        for epoch in range(epochs):
            # reset batch losses
            train_loss_batch = []
            for batch in loader_train:
                # output of autoencoder
                predicted = self.forward(batch)
                # calculating the loss
                batch_loss = loss_function(predicted, batch)
                # The gradients are set to zero,
                optimiser.zero_grad()
                # the the gradient is computed and stored.
                batch_loss.backward()
                # .step() performs parameter update
                optimiser.step()
                # Storing the losses for epoch
                train_loss_batch.append(float(batch_loss))

            validate_loss_batch = []
            with torch.set_grad_enabled(False):
                for batch in loader_validate:
                    predicted = self.forward(batch)
                    batch_loss = loss_function(predicted, batch)
                    validate_loss_batch.append(float(batch_loss))

            train_loss_epoch = sum(train_loss_batch)/len(train_loss_batch)
            validate_loss_epoch = sum(validate_loss_batch)/len(validate_loss_batch)
            losses.append([train_loss_epoch, ])
            print('Train: {:.4};  Validation: {:.4}'.format(train_loss_epoch, validate_loss_epoch))

            if epoch % 100 == 0:
                # Latent representation
                encoding = self.forward(loader_train, get_encoded=True)
                print(os.path.join(arg.out, "AE_epoch_" + str(epoch) + ".tsv"))
                # Save results
                np.savetxt(
                  os.path.join(arg.out, "AE_epoch_" + str(epoch) + ".tsv"),
                  encoding.cpu().detach().numpy(),
                  delimiter="\t"
                  )


if __name__ == "__main__":
    try:
        os.makedirs(arg.out)
    except FileExistsError:
        print("Directionary already exist for this training run, results will be overwritten")

    # Save input settings
    with open(arg.out + '/settings.json', 'w') as f:
        json.dump(arg.__dict__, f, indent=2)

    print("Loading features")
    # Load contig features
    features = pd.read_csv(arg.feature_file, sep='\t', header=0).values
    # features = pd.read_csv("/shared-nfs/SH/samples/zymo/mod_detection_test1/event_features.tsv", sep='\t', header=0).values
    features = features.astype("float")

    train, validate = torch.utils.data.random_split(features, [math.floor(features.shape[0]*(1-arg.val_percent)), math.ceil(features.shape[0]*(arg.val_percent))])
    train = train.dataset[1:(train.dataset.shape[0] - (train.dataset.shape[0] % arg.batch_size) + 1), ]
    validate = validate.dataset[1:(validate.dataset.shape[0] - (validate.dataset.shape[0] % arg.batch_size) + 1)]

    train = torch.from_numpy(train)
    train = train.float()
    train = torch.reshape(train, (-1, 2, arg.event_frame_size+1))

    validate = torch.from_numpy(validate)
    validate = validate.float()
    validate = torch.reshape(validate, (-1, 2, arg.event_frame_size+1))


    print("Initialising model")
    # Model Initialization
    model = AE_conv(
        input_shape=list(train.shape),
        hidden_nodes=arg.hidden_nodes,
        latent_nodes=arg.latent_nodes,
        kernel=arg.kernel_size,
        stride=1
    )

    print("Inititalising features")
    # Load dataset to train
    loader_train = torch.utils.data.DataLoader(
        dataset=train,
        batch_size=arg.batch_size,
        shuffle=False
    )
    loader_validate = torch.utils.data.DataLoader(
        dataset=validate,
        batch_size=arg.batch_size,
        shuffle=False
    )

    print("Training model")
    # model training
    model.train(
        loader_train,
        loader_validate,
        epochs=arg.epochs,
        learning_rate=arg.rate
    )

    # Latent representation
    encoding = model(loader_train, get_encoded=True)

    # Save results
    np.savetxt(arg.out + "/epoch_" + str(arg.epochs) + ".tsv", encoding.cpu().detach().numpy(), delimiter="\t")

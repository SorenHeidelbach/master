import argparse
# Parse arguments
# Instantiate the parser
parser = argparse.ArgumentParser(description='Embed multidimentional feature matrix into fewer latent representations with autoencoder')
# Required positional arguments
parser.add_argument('feature_file', metavar="features",
                    help='CSV file containing feature matrix')
parser.add_argument('out',
                    help='Directionary to save results too')
# Optional arguments
parser.add_argument('--latent_nodes', nargs=1, type=int, default=42, metavar='latent',
                    help="Number of latent nodes in bottleneck layer")
parser.add_argument('--hidden_nodes', nargs=1, type=int, default=200, metavar='hidden',
                    help="Number of nodes in the hidden layers of the encoder and decoder")
parser.add_argument('-t', '--threads', nargs=1, type=int, default=10,
                    help="Threads to allocate training of autoencoder")
parser.add_argument('--epochs', nargs=1, type=int, default=500,
                    help="Epochs to train autoencoder for")
parser.add_argument('--rate', nargs=1, type=int, default=0.01, metavar='rate',
                    help="Training rate")
parser.add_argument('--batch_size', nargs=1, type=int, default=100, metavar='batch',
                    help="Size of the batches in each epoch")
parser.add_argument('--seed', nargs=1, type=int, default=42,
                    help="Seed for pytorch and numpy")
# Load arguments
arg = parser.parse_args()

import torch
from torch import nn
torch.manual_seed(arg.seed)
torch.set_num_threads(arg.threads)
import pandas as pd
import numpy as np
np.random.seed(arg.seed)
import os, time, json, datetime


class AE(nn.Module):
    def __init__(self, input_shape, latent_nodes, hidden_nodes):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(in_features=input_shape, out_features=hidden_nodes),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_nodes),
            nn.Linear(in_features=hidden_nodes, out_features=hidden_nodes),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_nodes),
            nn.Linear(in_features=hidden_nodes, out_features=latent_nodes),
            nn.ReLU()
        )

        self.decoder = nn.Sequential(
            nn.Linear(in_features=latent_nodes, out_features=hidden_nodes),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_nodes),
            nn.Linear(in_features=hidden_nodes, out_features=hidden_nodes),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_nodes),
            nn.Linear(in_features=hidden_nodes, out_features=input_shape)

        )

    def forward(self, x, get_encoded=False):
        if get_encoded:
            x = self.encoder(x)
        else:
            x = self.encoder(x)
            x = self.decoder(x)
        return x

    def train(self, loader, epochs, learning_rate):
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
            batch_losses = []
            for batch in loader:
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
                batch_losses.append(float(batch_loss))

            epoch_loss = sum(batch_losses)/len(batch_losses)
            losses.append(epoch_loss)
            print("Epoch:" + str(epoch) + "   Loss:" + str(float(epoch_loss)))
            outputs.append((epoch, batch, predicted))

            if epoch % 100 == 0:
                # Latent representation
                encoding = self.forward(features, get_encoded=True)
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
    with open(arg.out + '/settings.txt', 'w') as f:
        json.dump(arg.__dict__, f, indent=2)

    print("Loading features")
    # Load contig features
    features = pd.read_csv(arg.feature_file, sep='\t').values
    features = features[:, 1:] # drop contig names
    features = features.astype("float")
    features = torch.from_numpy(features)
    features = features.float()

    print("Initialising model")
    # Model Initialization
    model = AE(
        input_shape=features.shape[1],
        hidden_nodes=arg.hidden_nodes,
        latent_nodes=arg.latent_nodes
    )

    print("Inititalising features")
    # Load dataset to train
    loader = torch.utils.data.DataLoader(
        dataset=features,
        batch_size=arg.batch_size,
        shuffle=False
    )

    print("Training model")
    # model training
    model.train(
        loader,
        epochs=arg.epochs,
        learning_rate=arg.rate
    )

    # Latent representation
    encoding = model(features, get_encoded=True)

    # Save results
    np.savetxt(arg.out + "/binnning_AE_representation.tsv", encoding.cpu().detach().numpy(), delimiter="\t")

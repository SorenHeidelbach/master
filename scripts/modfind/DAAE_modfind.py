#!/usr/bin/env python3

import argparse
# Parse arguments
# Instantiate the parser
parser = argparse.ArgumentParser(description='Learning categorical and style representation of current difference features from Nanodisco with Dual-AAE')
# Required positional arguments
parser.add_argument('feature_file', metavar="features",
                    help='tsv file containing abnomality feature matrix')
parser.add_argument('out',
                    help='Directionary to save results to')
# Optional arguments
parser.add_argument('--latent_nodes', type=int, default=2, metavar='latent',
                    help="Number of latent nodes in bottleneck style layer")
parser.add_argument('--categorical_nodes', type=int, default=50, metavar='latent',
                    help="Number of latent nodes in bottleneck categorical layer")

parser.add_argument('--hidden_nodes', type=int, default=100, metavar='hidden',
                    help="Number of nodes in the hidden layers of the encoder and decoder")
parser.add_argument('--kernel_size', type=int, default=3, metavar='kernel_size',
                    help="Size of kernel in convolution layer")
parser.add_argument('-t', '--threads', type=int, default=10,
                    help="Threads to allocate training of autoencoder")
parser.add_argument('--epochs', type=int, default=100,
                    help="Epochs to train autoencoder for")
parser.add_argument('--rate', type=int, default=0.0001, metavar='rate',
                    help="Training rate")
parser.add_argument('--batch_size', type=int, default=500, metavar='batch',
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
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.autograd import Variable
from torch.utils.data import DataLoader

import numpy as np
import pandas as pd
import os, json, csv, math

torch.manual_seed(arg.seed)
torch.set_num_threads(arg.threads)
np.random.seed(arg.seed)
##################################
# Define Networks
##################################
class Reshape(nn.Module):
    def __init__(self, *args):
        super().__init__()
        self.shape = args

    def forward(self, x):
        return x.view(self.shape)


class Q(nn.Module):
    ''' front end part of discriminator and encoder part of autoencoder'''

    def __init__(self, input_shape, categorical_nodes=10, latent_nodes=2, hidden_nodes=50, kernel=3, stride=1):
        super().__init__()
        self.size_conv=input_shape[2] - 3 * (kernel - 1)
        self.main = nn.Sequential(
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
            nn.Linear(self.size_conv*16, hidden_nodes),
            nn.BatchNorm1d(hidden_nodes),
            nn.LeakyReLU()
        )
        self.lin_cat = nn.Linear(hidden_nodes, categorical_nodes)
        self.lin_style = nn.Linear(hidden_nodes, latent_nodes)

    def forward(self, x):
        output = self.main(x)
        y = F.softmax(self.lin_cat(output), dim=0)
        z_gauss = self.lin_style(output)
        return z_gauss, y


class P(nn.Module):
    ''' Decoder '''

    def __init__(self, input_shape, latent_nodes=14, hidden_nodes=100, kernel=3, stride=1):
        super().__init__()
        self.size_conv=input_shape[2] - 3 * (kernel - 1)
        self.main = nn.Sequential(
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

    def forward(self, x):
        x = x.view(x.size(0), -1)
        output = self.main(x)
        return output


class D_gauss(nn.Module):
    ''' Discriminator for gaussian distrubution onto style embedding '''

    def __init__(self, latent_nodes=2):
        super().__init__()
        self.lin1 = nn.Linear(latent_nodes, 100)
        self.lin2 = nn.Linear(100, 100)
        self.lin3 = nn.Linear(100, 1)

    def forward(self, x):
        x = F.dropout(self.lin1(x), p=0.2, training=self.training)
        x = F.relu(x)
        x = F.dropout(self.lin2(x), p=0.2, training=self.training)
        x = F.relu(x)
        return self.lin3(x)


class Trainer:
    def __init__(self, Q, P, D_gauss, batch_size=100, epochs=500):
        assert Q.lin_style.out_features == D_gauss.lin1.in_features
        self.style_nodes = Q.lin_style.out_features
        self.cat_nodes = Q.lin_cat.out_features

        self.q = Q
        self.p = P
        self.d = D_gauss

        self.batch_size = batch_size
        self.epochs = epochs

    def sample_categorical(self, n_classes=10):
        '''
         Sample from a categorical distribution
         of size batch_size and # of classes n_classes
         return: torch.autograd.Variable with the sample
        '''
        idx = np.random.randint(0, n_classes, self.batch_size)
        cat = np.eye(n_classes)[idx].astype('float32')
        cat = torch.from_numpy(cat)
        return Variable(cat), idx

    def zerograd(self):
        self.q.zero_grad()
        self.p.zero_grad()
        self.d.zero_grad()

    def d_entropy1(self, y):
        y1 = torch.mean(y, 0)
        y2 = torch.sum(-y1 * torch.log(y1 + 1e-5))
        return y2

    def d_entropy2(self, y):
        y1 = -y * torch.log(y + 1e-5)
        y2 = torch.sum(y1) / self.batch_size
        return y2

    def train(self, loaded_data, val_data):
        ##################################
        # Load data and create Data loaders
        ##################################
        assert loaded_data.dataset.shape[0] % self.batch_size == 0
        assert val_data.dataset.shape[0] % self.batch_size == 0

        # loss function
        mse_dis = nn.MSELoss()
        criterionQ_dis = nn.NLLLoss()

        RE_solver = optim.Adam([{'params': self.q.parameters()}, {'params': self.p.parameters()}], lr=arg.rate)
        scheduler_RE = optim.lr_scheduler.MultiStepLR(RE_solver, milestones=[20,40], gamma=0.1)

        I_solver = optim.Adam([{'params': self.q.parameters()}, {'params': self.p.parameters()}], lr=arg.rate)
        scheduler_I = optim.lr_scheduler.MultiStepLR(I_solver, milestones=[20,40], gamma=0.1)

        Q_generator_solver = optim.Adam([{'params': self.q.parameters()}], lr=arg.rate)
        scheduler_Q = optim.lr_scheduler.MultiStepLR(Q_generator_solver, milestones=[20,40], gamma=0.1)

        D_loss_solver = optim.Adam([{'params': self.d.parameters()}], lr=arg.rate)
        scheduler_D = optim.lr_scheduler.MultiStepLR(D_loss_solver, milestones=[20,40], gamma=0.1)

        h_solver = optim.Adam([{'params': self.q.parameters()}], lr=arg.rate)
        scheduler_H = optim.lr_scheduler.MultiStepLR(h_solver, milestones=[20,40], gamma=0.2)

        losses=[["H", "D", "G", "Reconstruction", "I"]]
        losses_val=[["H", "D", "G", "Reconstruction", "I"]]
        # train network
        for epoch in range(self.epochs):
            print('\n', 'Epoch:', epoch, '\n')

            # Set the networks in train mode (apply dropout when needed)
            self.q.train()
            self.p.train()
            self.d.train()

            style_full = []
            cat_full = []

            for X in loaded_data:

                # Init gradients
                self.zerograd()

                # D loss
                z_real_gauss = Variable(torch.randn(self.batch_size, self.style_nodes))
                z_fake_gauss, ignore = self.q.forward(X)

                D_real_gauss = self.d(z_real_gauss)
                D_fake_gauss = self.d(z_fake_gauss.resize(self.batch_size, self.style_nodes))
                D_loss = -torch.mean(D_real_gauss) + torch.mean(D_fake_gauss)

                D_loss.backward()
                D_loss_solver.step()

                # Init gradients
                self.zerograd()

                # G loss
                z_fake_gauss, ignore = self.q(X)
                D_fake_gauss = self.d(z_fake_gauss.resize(self.batch_size, self.style_nodes))

                G_loss = -torch.mean(D_fake_gauss)
                G_loss.backward()
                Q_generator_solver.step()

                # Init gradients
                self.zerograd()

                # H loss
                ignore, y = self.q(X)
                h1_loss = self.d_entropy2(y)
                h2_loss = -self.d_entropy1(y)
                h_loss = 0.5 * (h1_loss + h2_loss)
                h_loss.backward()
                h_solver.step()

                # Init gradients
                self.zerograd()

                # recon_x
                z, y = self.q(X)
                X_re = self.p(torch.cat((z, y), 1).resize(self.batch_size, self.style_nodes + self.cat_nodes))

                recon_loss = mse_dis(X_re, X)
                recon_loss.backward()
                RE_solver.step()

                # Init gradients
                self.zerograd()

                # recon y and d
                y, idx = self.sample_categorical(n_classes=self.cat_nodes)
                z = Variable(torch.randn(self.batch_size, self.style_nodes))
                class_ = torch.LongTensor(idx)
                target = Variable(class_)

                X_sample = self.p(torch.cat((z, y), 1).resize(self.batch_size, self.style_nodes + self.cat_nodes, 1, 1))
                z_recon, y_recon = self.q(X_sample)

                I_loss = criterionQ_dis(torch.log(y_recon), target) + mse_dis(z_recon[:, 1:self.style_nodes], z[:, 1:self.style_nodes])

                I_loss.backward()
                I_solver.step()

                # Init gradients
                self.zerograd()

                if epoch % 20 == 0:
                    # Latent representation
                    style, categorical = self.q(X)
                    style_full.append(style.detach().cpu().numpy())
                    cat_full.append(categorical.detach().cpu().numpy())
            with torch.set_grad_enabled(False):
                for X in val_data:
                    # Init gradients
                    self.zerograd()

                    # D loss
                    z_real_gauss = Variable(torch.randn(self.batch_size, self.style_nodes))
                    z_fake_gauss, ignore = self.q.forward(X)

                    D_real_gauss = self.d(z_real_gauss)
                    D_fake_gauss = self.d(z_fake_gauss.resize(self.batch_size, self.style_nodes))
                    D_loss_val = -torch.mean(D_real_gauss) + torch.mean(D_fake_gauss)

                    # Init gradients
                    self.zerograd()

                    # G loss
                    z_fake_gauss, ignore = self.q(X)
                    D_fake_gauss = self.d(z_fake_gauss.resize(self.batch_size, self.style_nodes))

                    G_loss_val = -torch.mean(D_fake_gauss)

                    # Init gradients
                    self.zerograd()

                    # H loss
                    ignore, y = self.q(X)
                    h1_loss = self.d_entropy2(y)
                    h2_loss = -self.d_entropy1(y)

                    h_loss_val = 0.5 * (h1_loss + h2_loss)

                    # Init gradients
                    self.zerograd()

                    # recon_x
                    z, y = self.q(X)
                    X_re = self.p(torch.cat((z, y), 1).resize(self.batch_size, self.style_nodes + self.cat_nodes))

                    recon_loss_val = mse_dis(X_re, X)

                    # Init gradients
                    self.zerograd()

                    # recon y and d
                    y, idx = self.sample_categorical(n_classes=self.cat_nodes)
                    z = Variable(torch.randn(self.batch_size, self.style_nodes))
                    class_ = torch.LongTensor(idx)
                    target = Variable(class_)

                    X_sample = self.p(
                        torch.cat((z, y), 1).resize(self.batch_size, self.style_nodes + self.cat_nodes, 1, 1))
                    z_recon, y_recon = self.q(X_sample)

                    I_loss_val = criterionQ_dis(torch.log(y_recon), target) + mse_dis(z_recon[:, 1:4], z[:, 1:4])


            # Test model
            self.q.eval()

            # Report training losses
            print('Train:      H: {:.4}; D: {:.4}; G: {:.4}; recon: {:.4}; I: {:.4}'.format(
                h_loss.item(), D_loss.item(), G_loss.item(), recon_loss.item(), I_loss.item()))
            losses.append([h_loss.item(), D_loss.item(), G_loss.item(), recon_loss.item(), I_loss.item()])

            # Report validation losses
            print('Validation: H: {:.4}; D: {:.4}; G: {:.4}; recon: {:.4}; I: {:.4}'.format(
                h_loss_val.item(), D_loss_val.item(), G_loss_val.item(), recon_loss_val.item(), I_loss_val.item()))
            losses_val.append([h_loss_val.item(), D_loss_val.item(), G_loss_val.item(), recon_loss_val.item(), I_loss_val.item()])


            if epoch % 20 == 0:
                # Save results
                # Style embedding
                np.savetxt(
                    os.path.join(arg.out + "/embedding_style_epoch_" + str(epoch) + ".tsv"),
                    np.concatenate(style_full),
                    delimiter="\t"
                  )
                # Categerical embedding
                np.savetxt(
                    os.path.join(arg.out + "/embedding_categorical_epoch_" + str(epoch) + ".tsv"),
                    np.concatenate(cat_full),
                    delimiter="\t"
                )
                # Training losses
                with open(os.path.join(arg.out + "/loss_training_epoch_" + str(epoch) + ".tsv"), 'w', newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerows(losses)
                # Validation losses
                with open(os.path.join(arg.out + "/loss_validation_epoch_" + str(epoch) + ".tsv"), 'w', newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerows(losses_val)
            # Update learning rate
            scheduler_RE.step()
            scheduler_I.step()
            scheduler_Q.step()
            scheduler_D.step()
            scheduler_H.step()
        return losses, losses_val

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
    Q_mod = Q(
        list(train.shape),
        latent_nodes=arg.latent_nodes,
        categorical_nodes=arg.categorical_nodes,
        hidden_nodes=arg.hidden_nodes,
        kernel=arg.kernel_size,
        stride=1
    )
    P_mod = P(
        list(train.shape),
        latent_nodes=arg.categorical_nodes + arg.latent_nodes
    )
    D_mod = D_gauss(
        latent_nodes=arg.latent_nodes
    )

    model = Trainer(
        Q=Q_mod,
        P=P_mod,
        D_gauss=D_mod,
        batch_size=arg.batch_size,
        epochs=arg.epochs
    )

    print("Inititalising features")
    # Load dataset to train
    loader = torch.utils.data.DataLoader(
        dataset=train,
        batch_size=arg.batch_size,
        shuffle=False
    )
    loader_val = torch.utils.data.DataLoader(
        dataset=validate,
        batch_size=arg.batch_size,
        shuffle=False
    )


    print("Training model")
    # model training
    loss, val_loss = model.train(loaded_data=loader, val_data=loader_val)

    # Save results
    style, cat = model.q(train)
    np.savetxt(
      os.path.join(arg.out + "/embedding_categorical_epoch_" + str(arg.epochs) + ".tsv"),
      cat.detach().cpu().numpy(),
      delimiter="\t"
    )
    np.savetxt(
      os.path.join(arg.out + "/embedding_style_epoch_" + str(arg.epochs) + ".tsv"),
      cat.detach().cpu().numpy(),
      delimiter="\t"
    )
    # Training losses
    with open(os.path.join(arg.out + "/loss_training_epoch_" + str(arg.epochs) + ".tsv"), 'w', newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(loss)
    # Validation losses
    with open(os.path.join(arg.out + "/loss_validation_epoch_" + str(arg.epochs) + ".tsv"), 'w', newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(val_loss)

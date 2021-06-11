import torch
import torchvision
import torchvision.transforms as transforms
import numpy as np
import math
import torch.nn as nn
import torch.nn.functional as F

device = torch.device("cpu")

def get_counts(countfile):
    npzFile = np.load(countfile)
    matrix = npzFile['count']
    npzFile.close()
    return matrix

def numeric_rep(frontierfile,countfile,K,k):
    fp = open(frontierfile)
    data = fp.read().splitlines()
    fp.close()
    n = len(data)
    bd = {'A':0, 'C':1, 'G':2, 'T':3, 'N': 0}
    x = []
    counts = get_counts(countfile)
    seq_list = []
    for j,seq in enumerate(data):
        seq_list.append(seq)
        x_part = np.zeros((4,K-k+1))
        for i,b in enumerate(seq[0:K-k+1]):
            bb = seq[k//2+i]
            x_part[bd[bb]][i] = math.log2(int(counts[j][i]))
            #x_part[bd[b]][i] = 1
        x.append(x_part)
    fp.close()
    X_test = np.array(x)
    X_test = torch.from_numpy(X_test).float()
    return X_test, seq_list
	
class frontierNet(nn.Module):
    def __init__(self, D_out):
        super(frontierNet, self).__init__()
        self.conv1 = nn.Conv1d(4, 9, 15, stride=1)
        self.pool = nn.MaxPool1d(4, stride=2)
        self.conv2 = nn.Conv1d(9, 12, 5, stride=1)
        self.fc1 = nn.Linear(12*6, 12*6)
        self.fc2 = nn.Linear(12*6, D_out)
        self.fc3 = nn.Linear(12*6, D_out)
        #self.prob = nn.Softmax(dim=1)
        
    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        #print (x.shape)
        x = x.view(-1, 12 * 6)
        #print (x.shape)
        x = self.fc1(x)
        x = self.fc2(x)
        #x = self.fc3(x)
        #x = self.prob(x)
        return x

#net = frontierNet(D_out)

def run(frontierfile,countfile,PATH,truefrontier):
	K = 75
	k = 21
	X_test, seq_list = numeric_rep(frontierfile,countfile,K,k)
	N = 64
	D_out = 4

	net = torch.load(PATH)
	net.eval()

	Y_pred = torch.argmax(net(X_test), dim=1)

	fp = open(truefrontier, "w")
	for i,c in enumerate(Y_pred):
		if c == 0:
			fp.write("%s\n" %seq_list[i])
	fp.close()
	
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('frontierfile', help='List of frontier candidates')
    parser.add_argument('countfile', help='Desired Path and Name of the countfile')
	parser.add_argument('trainedmodel', help='Path to the trained model for finding true frontier')
	parser.add_argument('truefrontier', help='Desired Path to the output that keeps only the true frontier sequences')
    
    args = parser.parse_args()
    
    frontierfile = args.frontierfile
	countfile = args.countfile
	PATH = args.trainedmodel
	truefrontier = args.truefrontier
    run(frontierfile, countfile, PATH, truefrontier)
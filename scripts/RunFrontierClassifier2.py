import torch
import torchvision
import torchvision.transforms as transforms
import numpy as np
import math
import torch.nn as nn
import torch.nn.functional as F

device = torch.device("cpu")

def numeric_rep(truefrontier, K,k):
    fp = open(truefrontieru)
    data = fp.read().splitlines()
    fp.close()
    n = len(data)
    bd = {'A':0, 'C':1, 'G':2, 'T':3, 'N': 0}
    x, seq_list = [], []
    cls_dict = {"LINE":0, "SINE_Alu":1, "SINE_B2":2, "LTR":3}
    for j,seq in enumerate(data):
        c = cls_dict[fam]
        seq_list.append(seq)
        x_part = np.zeros((4,K))
        for i,b in enumerate(seq[0:K]):
            x_part[bd[b]][i] = 1
        x.append(x_part)
        y.append(c)
    fp.close()
    X_test = np.array(x)
    Y_train = np.array(y)
    X_test = torch.from_numpy(X_train).float()
    return X_test, seq_list
	
class frontierNet(nn.Module):
    def __init__(self, D_out):
        super(frontierNet, self).__init__()
        self.conv1 = nn.Conv1d(4, 9, 15, stride=1)
        self.pool = nn.MaxPool1d(4, stride=2)
        self.conv2 = nn.Conv1d(9, 12, 5, stride=1)
        self.fc1 = nn.Linear(12*11, 12*11)
        self.fc2 = nn.Linear(12*11, D_out)
        self.fc3 = nn.Linear(12*11, D_out)
        #self.prob = nn.Softmax(dim=1)
        
    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        #print (x.shape)
        x = x.view(-1, 12 * 11)
        #print (x.shape)
        x = self.fc1(x)
        x = self.fc2(x)
        #x = self.fc3(x)
        #x = self.prob(x)
        return x
		
def run(PATH,truefrontier,truefrontier_classified):
	K = 75
	k = 21
	N = 64
	D_out = 4
	X_test, seq_list = numeric_rep(truefrontier,K,k)

	net = torch.load(PATH)
	net.eval()

	cls_list = ["LINE", "SINE_Alu", "SINE_B2", "LTR"]
	Y_pred = torch.argmax(net(X_test), dim=1)

	fp = open(truefrontier_classified, "w")
	for i,c in enumerate(Y_pred):
		fp.write("%s,%s\n" %(seq_list[i],cls_list[c]))
	fp.close()
	
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
	parser.add_argument('truefrontier', help='Path to the true frontiers - output of 1st classifier')
	parser.add_argument('trainedmodel', help='Path to the trained model for classifier 2')
    parser.add_argument('truefrontier_classified', help='Desired Path to the true frontiers with TE type classified')
    args = parser.parse_args()
    
    truefrontier = args.truefrontier
	PATH = args.trainedmodel
	truefrontier_classified = args.truefrontier_classified
    run(truefrontier, PATH, truefrontier_classified)